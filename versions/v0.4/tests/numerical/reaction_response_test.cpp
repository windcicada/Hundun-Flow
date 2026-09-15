// SPDX-License-Identifier: Apache-2.0
#include "../../src/core_response.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

using namespace hundun::v04;
int main() {
  detail::IntervalReactionResponse response;
  bool pass = bool(response.prepare(2, 3, 2, 1e-11, 1e-20));
  const auto bytes = response.bytes();
  std::array<double,3> initial{.3,.2,.5}, first{.2,.3,.5};
  detail::IntervalReactionResponse::Report report;
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= !report.reused;
  const auto original = first;
  // A mass- and element-preserving response perturbation is below the
  // specified coupling budget. Retain the same whole reaction increment.
  first[0] += 1e-14; first[1] -= 1e-14;
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= report.reused && report.error_ratio > 0 && report.error_ratio <= 1;
  pass &= first == original && response.bytes() == bytes;
  // The same endpoint perturbation exceeds a short physical interval's
  // rate-based budget. Resolve it instead of reusing the longer-step value.
  first[0] += 1e-14; first[1] -= 1e-14;
  const auto refined = first;
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report,1e-5));
  pass &= !report.reused && first == refined;
  for (double scale : {-1., 2., std::numeric_limits<double>::quiet_NaN()}) {
    pass &= !response.select(0,{initial.data(),3},{first.data(),3},report,scale);
    pass &= first == refined;
  }
  response.reset();
  first = original;
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  // Transport changes the local composition. Replay its reaction increment
  // on the new transport endpoint instead of freezing the final composition.
  initial = {.29,.21,.5}; first = {.19+1e-14,.31-1e-14,.5};
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= report.reused && std::abs(first[0]-.19)<1e-16 &&
          std::abs(first[1]-.31)<1e-16;
  // Resolve a physical response change immediately and atomically.
  first = {.18,.32,.5};
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= !report.reused && first[0] == .18 && first[1] == .32;
  // Each cell owns its response. Reset discards every failed-attempt response.
  first = {.18+1e-14,.32-1e-14,.5};
  pass &= bool(response.select(1,{initial.data(),3},{first.data(),3},report));
  pass &= !report.reused;
  response.reset();
  pass &= bool(response.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= !report.reused;
  // A weak component uses its absolute error budget and stays admissible.
  detail::IntervalReactionResponse trace;
  pass &= bool(trace.prepare(1,3,2,1e-11,1e-20));
  initial = {1e-21,.2,.8}; first = {0.,.2,.8};
  pass &= bool(trace.select(0,{initial.data(),3},{first.data(),3},report));
  initial[0]=0.;
  pass &= bool(trace.select(0,{initial.data(),3},{first.data(),3},report));
  pass &= !report.reused && first[0]==0.;
  // Invalid fresh outputs leave both the proposed endpoint and cache intact.
  const auto saved=first;
  first[0]=std::numeric_limits<double>::quiet_NaN();
  pass &= !trace.select(0,{initial.data(),3},{first.data(),3},report);
  pass &= std::isnan(first[0]) && first[1]==saved[1];
  pass &= !trace.prepare(SIZE_MAX,3,2,1e-11,1e-20);
  if(!pass) std::cerr << "interval response coupling budget/reset/admissibility failure\n";
  return pass ? 0 : 1;
}

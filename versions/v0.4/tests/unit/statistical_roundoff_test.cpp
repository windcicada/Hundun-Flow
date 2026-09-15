// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_statistical_detail.hpp"
#include <cstring>
#include <iostream>

using hundun::v04::detail::close_statistical_composition;
int main() {
  unsigned failures=0;
  const auto check=[&](bool ok){failures+=!ok;};
  const double eps=std::numeric_limits<double>::epsilon();
  for (std::size_t dependent=0;dependent<3;++dependent) {
    const auto a=(dependent+1)%3,b=(dependent+2)%3;
    std::array<double,3> y{};
    y[a]=-3.4538625640398644e-27;y[b]=.6;
    check(close_statistical_composition({y.data(),y.size()},dependent));
    check(y[a]==0 && y[b]==.6 && y[dependent]==.4);
    y[a]=.5+eps;y[b]=.5;
    check(close_statistical_composition({y.data(),y.size()},dependent));
    check(y[dependent]>=0 && static_cast<long double>(y[a])+y[b]<=1);
    check(std::abs(y[a]-(.5+eps))+std::abs(y[b]-.5)<=2*eps);
    y[a]=.2;y[b]=.3;
    check(close_statistical_composition({y.data(),y.size()},dependent));
    check(y[a]==.2 && y[b]==.3);
    for(double invalid:{-1e-12,1+1e-12,std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}) {
      y[a]=invalid;y[b]=.2;y[dependent]=.7;
      const auto before=y;
      check(!close_statistical_composition({y.data(),y.size()},dependent));
      check(std::memcmp(before.data(),y.data(),sizeof(y))==0);
    }
    y[a]=.6;y[b]=.4+1e-12;
    const auto before=y;
    check(!close_statistical_composition({y.data(),y.size()},dependent));
    check(std::memcmp(before.data(),y.data(),sizeof(y))==0);
  }
  // Many individually small negatives must obey one total repair budget.
  std::array<double,40> many{};many.fill(-eps);
  const auto before=many;
  check(!close_statistical_composition({many.data(),many.size()},39));
  check(std::memcmp(before.data(),many.data(),sizeof(many))==0);
  for(unsigned units=1;units<=40;++units) {
    std::array<double,3> y{-double(units)*eps,.123456789,0};
    const auto before=y;
    const long double dependent=1-static_cast<long double>(y[0])-y[1];
    if(close_statistical_composition({y.data(),y.size()},2))
      check(std::abs(static_cast<long double>(y[0])-before[0])+
            std::abs(static_cast<long double>(y[1])-before[1])+
            std::abs(static_cast<long double>(y[2])-dependent)<=64*eps);
    else check(std::memcmp(before.data(),y.data(),sizeof(y))==0);
  }
  std::cerr<<"statistical roundoff failures="<<failures<<'\n';
  return failures?1:0;
}

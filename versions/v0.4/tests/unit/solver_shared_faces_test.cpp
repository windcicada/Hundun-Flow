// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "../../src/solver_shared_faces_detail.hpp"
#include <cmath>
#include <iostream>
#include <vector>
using namespace hundun::v04;
int main() {
  constexpr Int3 cells{19,13,11};
  constexpr KernelBox box{{1,2,1},{17,10,9}};
  std::vector<std::uint8_t> active(cells.x*cells.y*cells.z,1U);
  double offset=0.0;
  const auto value=[&](CartesianAxis axis,Int3 f) {
    return offset+std::sin(0.17*f.x+0.23*f.y+0.37*f.z+static_cast<unsigned>(axis));
  };
  bool passed=true;
  for(bool masked:{false,true}) {
    if(masked) for(std::size_t i=0;i<active.size();++i) active[i]=(i%7!=0 && i%13!=0);
    offset+=0.1;
    std::vector<unsigned> visits(active.size());
    std::size_t calls=0, visited=0;
    auto status=detail::shared_cell_faces(cells,box,{active.data(),active.size()},
      [&](CartesianAxis a,Int3 f,double& out) { ++calls;out=value(a,f);return Status{}; },
      [&](Int3 c,const std::array<double,6U>& faces) {
        ++visited;
        passed &= ++visits[(c.z*cells.y+c.y)*cells.x+c.x]==1U;
        for(unsigned a=0;a<3;++a) for(unsigned side=0;side<2;++side) {
          auto f=c;(a==0 ? f.x : a==1 ? f.y : f.z)+=side;
          passed &= faces[2*a+side]==value(static_cast<CartesianAxis>(a),f);
        }
        return Status{};
      });
    const auto expected=[&]() { std::size_t n=0;for(int z=1;z<10;++z) for(int y=2;y<12;++y) for(int x=1;x<18;++x) n+=active[(z*cells.y+y)*cells.x+x]!=0;return n; }();
    passed &= status && visited==expected && calls<5*visited;
    std::cout<<"shared_faces masked="<<masked<<" cells="<<visited<<" evaluations="<<calls<<" naive="<<6*visited<<'\n';
  }
  unsigned consumed=0;
  auto status=detail::shared_cell_faces(cells,box,{},
    [](CartesianAxis,Int3,double&) {return Status{StatusCode::numerical_failure,456U};},
    [&](Int3,const std::array<double,6U>&) {++consumed;return Status{};});
  passed &= !status && status.detail==456U && consumed==0U;
  return passed ? 0 : 1;
}

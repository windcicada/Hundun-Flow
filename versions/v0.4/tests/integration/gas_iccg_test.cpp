// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_iccg.hpp"
#include <iostream>
using namespace hundun::v04;
namespace {
struct Field {
  std::vector<double> data;
  FieldView view;
  Field(FieldId id,Int3 n):data(std::size_t(n.x)*n.y*n.z) {
    view.base=data.data();view.interior=n;view.components=1;
    view.stride_y=n.x;view.stride_z=std::size_t(n.x)*n.y;
    view.component_stride=data.size();view.field=id;view.revision=1;
    view.storage_identity=100+id;view.revision_domain=99;
  }
};
}
int main() {
  constexpr Int3 cells{3,2,2};
  std::vector<detail::ColdPressureRow> rows(12);
  for(std::size_t i=0;i<rows.size();++i) {
    rows[i].diagonal=20+i*.25;
    // The preconditioner owns an explicit L D L^T even for a nonsymmetric
    // input. Its certificate proves M SPD; PCG must also prove A SPD.
    rows[i].neighbour={.25+i/32.,3.,.5+i/32.,4.,.75+i/32.,5.};
  }
  LinearIdentity identity{11,12,13,14,15};
  detail::ColdPressureIccg pc(rows,cells,identity);
  Field x(1,cells),y(2,cells),mx(3,cells),my(4,cells);
  bool pass=!pc.apply(as_const(x.view),mx.view,0);
  auto status=pc.prepare();pass=pass && status;
  pass=pass && pc.certificate().preconditioner_class==LinearPreconditionerClass::fixed_spd;
  for(std::size_t i=0;i<rows.size();++i) {x.data[i]=(int(i)-5)/8.;y.data[i]=(int(i)%3-1)/4.;}
  status=pc.apply(as_const(x.view),mx.view,0);pass=pass && status;
  status=pc.apply(as_const(y.view),my.view,1);pass=pass && status;
  double xmy=0,ymx=0,xmx=0;
  for(std::size_t i=0;i<rows.size();++i) {xmy+=x.data[i]*my.data[i];ymx+=y.data[i]*mx.data[i];xmx+=x.data[i]*mx.data[i];}
  pass=pass && std::abs(xmy-ymx)<1e-15 && xmx>0;
  // Reconstruct the dense M=(D-L) D^-1 (D-L)^T independently and check M z=b.
  const auto pivots=pc.inverse_pivots();
  std::array<std::array<double,12>,12> lower{};
  for(std::size_t i=0;i<12;++i) {
    lower[i][i]=1/pivots.data[i];
    if(i%3) lower[i][i-1]=-rows[i].neighbour[0];
    if((i/3)%2) lower[i][i-3]=-rows[i].neighbour[2];
    if(i/6) lower[i][i-6]=-rows[i].neighbour[4];
  }
  double error=0;
  for(std::size_t i=0;i<12;++i) {
    double sum=0;
    for(std::size_t j=0;j<12;++j) {
      double m=0;for(std::size_t k=0;k<12;++k)m+=lower[i][k]*pivots.data[k]*lower[j][k];
      sum+=m*mx.data[j];
    }
    error=std::max(error,std::abs(sum-x.data[i]));
  }
  pass=pass && error<1e-14;
  rows[1].diagonal=0;
  const auto before=mx.data;
  pass=pass && !pc.prepare() && pc.inverse_pivots().size==0 &&
      !pc.apply(as_const(x.view),mx.view,0) && mx.data==before;
  detail::ColdPressureRow original{{1,2,3,4,5,6},21,7},scaled;
  status=detail::scale_iccg_row(original,.125,scaled);
  pass=pass && status && scaled.diagonal==21*.125 && scaled.rhs==7*.125;
  for(unsigned f=0;f<6;++f)pass=pass && scaled.neighbour[f]==original.neighbour[f]*.125;
  for(double volume : {0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::max()}) {
    const auto before=scaled;
    pass=pass && !detail::scale_iccg_row(original,volume,scaled) &&
        before.diagonal==scaled.diagonal && before.neighbour==scaled.neighbour && before.rhs==scaled.rhs;
  }
  std::cout<<"ICCG transpose defect="<<std::abs(xmy-ymx)<<" dense factor residual="<<error<<" pass="<<pass<<'\n';
  return pass ? 0:1;
}

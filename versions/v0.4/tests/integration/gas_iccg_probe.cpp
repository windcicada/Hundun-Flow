// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_iccg.hpp"
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace hundun::v04;
namespace {
struct Field {
  std::vector<double> data;
  FieldView view;
  Field(FieldId id,Int3 n,unsigned components=1)
      :data(std::size_t(n.x)*n.y*n.z*components) {
    view.base=data.data();view.interior=n;view.components=components;
    view.stride_y=n.x;view.stride_z=std::size_t(n.x)*n.y;
    view.component_stride=std::size_t(n.x)*n.y*n.z;
    view.field=id;view.revision=1;view.storage_identity=100+id;view.revision_domain=99;
  }
};
class Matrix final : public LinearOperator {
 public:
  const std::vector<detail::ColdPressureRow>& rows;
  Int3 n;std::array<int,3> periodic;
  LinearIdentity identity;
  Matrix(const std::vector<detail::ColdPressureRow>& r,Int3 cells,
      std::array<int,3> p,LinearIdentity id):rows(r),n(cells),periodic(p),identity(id) {}
  Int3 neighbour(Int3 c,unsigned f) const {
    auto& v=f/2==0 ? c.x : f/2==1 ? c.y : c.z;
    const auto extent=f/2==0 ? n.x : f/2==1 ? n.y : n.z;
    v+=f%2 ? 1 : -1;
    if(periodic[f/2])v=(v+extent)%extent;
    return c;
  }
  std::size_t flat(Int3 c) const {return std::size_t(c.x)+n.x*(c.y+std::size_t(n.y)*c.z);}
  // This whole-grid probe owns all rows. Exact binary coefficients make the
  // graph/M-matrix proof an independent gate for its SPD operator certificate.
  bool admissible() const {
    std::vector<std::size_t> parent(rows.size());std::iota(parent.begin(),parent.end(),0);
    std::vector<bool> ground(rows.size());
    const auto root=[&](std::size_t i) {while(parent[i]!=i)i=parent[i];return i;};
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      Int3 c{x,y,z};auto i=flat(c);long double sum=0;
      if(!(rows[i].diagonal>0) || !std::isfinite(rows[i].diagonal))return false;
      for(unsigned f=0;f<6;++f) {
        const double a=rows[i].neighbour[f];if(!std::isfinite(a)||a<0)return false;
        sum+=a;if(a==0)continue;
        const auto nb=neighbour(c,f);
        if(nb.x<0 || nb.x>=n.x || nb.y<0 || nb.y>=n.y || nb.z<0 || nb.z>=n.z)return false;
        const auto j=flat(nb);
        if(a!=rows[j].neighbour[f^1])return false;
        parent[root(j)]=root(i);
      }
      if(rows[i].diagonal<sum)return false;
      ground[i]=rows[i].diagonal>sum;
    }
    for(std::size_t i=0;i<rows.size();++i) if(ground[i])ground[root(i)]=true;
    for(std::size_t i=0;i<rows.size();++i)if(!ground[root(i)])return false;
    return true;
  }
  LinearOperatorCertificate certificate() const noexcept override {
    return {identity,19,n,LinearOperatorClass::spd};
  }
  Status apply(FieldView in,FieldView out) const noexcept override {
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      const Int3 c{x,y,z};const auto& row=rows[flat(c)];
      double v=row.diagonal*in.unchecked(c,0);
      for(unsigned f=0;f<6;++f)if(row.neighbour[f]!=0)v-=row.neighbour[f]*in.unchecked(neighbour(c,f),0);
      out.unchecked(c,0)=v;
    }
    return {};
  }
};
bool run(Int3 n,std::array<int,3> periodic) {
  if(n.x<1 || n.y<1 || n.z<1 || n.x>64 || n.y>64 || n.z>64)return false;
  const auto count=std::size_t(n.x)*n.y*n.z;
  std::vector<detail::ColdPressureRow> original(count),rows(count);
  std::vector<double> volumes(count);
  for(std::size_t i=0;i<count;++i) {
    auto& row=original[i];std::cin>>volumes[i]>>row.diagonal;
    for(auto& a:row.neighbour)std::cin>>a;
    std::cin>>row.rhs;
    if(!std::cin || !detail::scale_iccg_row(row,volumes[i],rows[i]))return false;
  }
  LinearWorkspaceRequirements req;
  auto status=make_linear_workspace_requirements(LinearAlgorithm::pcg,n,0,0,
      ReductionMode::mpi_allreduce,17,req);
  Field b(1,n),x(2,n),r(3,n),vectors(4,n,req.vector_slots),scalars(5,{int(req.scalar_doubles),1,1});
  scalars.view.storage_identity=vectors.view.storage_identity;
  SolverWorkspace work;
  if(status)status=SolverWorkspace::bind(req,vectors.view,scalars.view,work);
  ReductionEngine reductions;
  if(status)status=ReductionEngine::compile(MPI_COMM_SELF,ReductionMode::mpi_allreduce,req.reduction_capacity,reductions);
  LinearIdentity identity{11,12,13,work.fingerprint(),15};
  Matrix op(rows,n,periodic,identity);
  if(!status || !op.admissible())return false;
  detail::ColdPressureIccg pc(rows,n,identity);
  status=pc.prepare();if(!status)return false;
  for(std::size_t i=0;i<count;++i)b.data[i]=rows[i].rhs;
  const LinearSolveControl control{1e-12,1e-13,2000,64,0};
  const auto solved=solve_pcg(op,pc,{as_const(b.view),x.view,identity,control},work,reductions);
  if(!solved.status)return false;
  status=op.apply(x.view,r.view);if(!status)return false;
  double residual=0;
  for(std::size_t i=0;i<count;++i)residual=std::max(residual,std::abs((b.data[i]-r.data[i])/volumes[i]));
  std::cout<<std::setprecision(17)<<solved.iterations<<' '<<residual<<'\n';
  const auto inv=pc.inverse_pivots();
  for(std::size_t i=0;i<count;++i) {
    std::cout<<x.data[i]<<' '<<inv.data[i]<<' '<<rows[i].diagonal;
    for(auto a:rows[i].neighbour)std::cout<<' '<<a;
    std::cout<<' '<<rows[i].rhs<<'\n';
  }
  return true;
}
}
int main(int argc,char** argv) {
  if(MPI_Init(&argc,&argv)!=MPI_SUCCESS)return 2;
  bool okay=true;Int3 n;std::array<int,3> periodic{};
  while(std::cin>>n.x>>n.y>>n.z>>periodic[0]>>periodic[1]>>periodic[2])
    if(!run(n,periodic)){okay=false;break;}
  MPI_Finalize();return okay ? 0:1;
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "../../src/solver_cold.hpp"
#include <array>
#include <iostream>
#include <vector>
using namespace hundun::v04;
namespace {
struct Field {
  std::vector<double> data;FieldView view;
  Field(FieldId id,int cells,unsigned components=1):data(std::size_t(cells)*components) {
    view.base=data.data();view.interior={cells,1,1};view.components=components;
    view.stride_y=view.stride_z=cells;view.component_stride=cells;
    view.field=id;view.revision=1;view.storage_identity=1000+id;view.revision_domain=91;
  }
};
class Identity final : public LinearOperator {
 public:
  explicit Identity(LinearIdentity identity):identity_(identity) {}
  LinearOperatorCertificate certificate() const noexcept override {return {identity_,11,{3,1,1},LinearOperatorClass::nonsymmetric};}
  Status apply(FieldView x,FieldView y) const noexcept override {
    for(int i=0;i<3;++i)y.unchecked({i,0,0},0)=x.unchecked({i,0,0},0);
    return {};
  }
 private:LinearIdentity identity_;
};
class UnitPreconditioner final : public LinearPreconditioner {
 public:
  explicit UnitPreconditioner(LinearIdentity id):id_(id) {}
  LinearPreconditionerCertificate certificate() const noexcept override {return {id_,12,LinearPreconditionerClass::fixed_general};}
  Status apply(ConstFieldView x,FieldView y,std::uint32_t) noexcept override {
    for(int i=0;i<3;++i)y.unchecked({i,0,0},0)=x.unchecked({i,0,0},0);return {};
  }
 private:LinearIdentity id_;
};
bool run(LinearAlgorithm algorithm,int rank,int size) {
  const unsigned restart=algorithm==LinearAlgorithm::fgmres ? 8 : 0;
  LinearWorkspaceRequirements req;
  auto status=make_linear_workspace_requirements(algorithm,{3,1,1},0,restart,ReductionMode::mpi_allreduce,71,req);
  if(!status){std::cerr<<"audit_setup rank="<<rank<<" status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  Field b(1,3),x(2,3),rho(3,3),p(4,3),residual(5,3),vectors(6,3,req.vector_slots),scalars(7,req.scalar_doubles);
  scalars.view.storage_identity=vectors.view.storage_identity;
  SolverWorkspace workspace;status=SolverWorkspace::bind(req,vectors.view,scalars.view,workspace);
  ReductionEngine reductions;
  if(status)status=ReductionEngine::compile(MPI_COMM_WORLD,ReductionMode::mpi_allreduce,req.reduction_capacity,reductions);
  if(!status){std::cerr<<"audit_setup rank="<<rank<<" status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  LinearIdentity identity{11,12,13,workspace.fingerprint(),15};
  Identity op(identity);UnitPreconditioner pc(identity);
  for(auto& v:rho.data)v=1;
  const auto initialize=[&] {
    b.data[0]=x.data[0]=1e6;b.data[1]=1;x.data[1]=rank==size-1 ? 1.-1e-8 : 1.;
    b.data[2]=x.data[2]=0;rho.data[1]=1e-3;
  };
  initialize();
  const LinearSolveControl control{1e-13,1e-13,100,4,restart};
  const auto solve=[&](LinearConvergenceAudit* audit) {
    LinearSolveInvocation call{as_const(b.view),x.view,identity,control,audit};
    return algorithm==LinearAlgorithm::fgmres ? solve_fgmres(op,pc,call,workspace,reductions)
        : solve_bicgstab(op,pc,call,workspace,reductions);
  };
  const auto baseline=solve(nullptr);
  detail::ColdPressureContinuityAudit audit(as_const(rho.view),as_const(p.view),1e5,1.,{},8123);
  for(int i=0;i<3;++i)residual.data[i]=b.data[i]-x.data[i];
  LinearConvergenceAuditResult before;
  status=audit.evaluate(as_const(x.view),as_const(residual.view),reductions,before);
  bool okay=status && baseline.status && baseline.iterations==0 && !before.accepted && before.metric>1e-6;
  initialize();const auto corrected=solve(&audit);
  for(int i=0;i<3;++i)residual.data[i]=b.data[i]-x.data[i];
  LinearConvergenceAuditResult after;
  status=audit.evaluate(as_const(x.view),as_const(residual.view),reductions,after);
  okay=okay && status && corrected.status && corrected.iterations>0 && corrected.convergence_rejections>=1 &&
      after.accepted && x.data==b.data && corrected.true_residual_limit==baseline.true_residual_limit &&
      corrected.relative_tolerance==control.relative_tolerance && corrected.absolute_tolerance==control.absolute_tolerance;
  // Solid rows remain in the canonical linear solve; the fluid-scale audit
  // excludes their thermophysical state. Only the final rank carries a defect.
  std::array<std::uint8_t,3> active{1,0,1};
  detail::ColdPressureContinuityAudit masked(as_const(rho.view),as_const(p.view),1e5,1.,{active.data(),3},8124);
  if(rank==size-1) {rho.data[1]=std::numeric_limits<double>::quiet_NaN();residual.data[1]=1;}
  LinearConvergenceAuditResult ignored;
  status=masked.evaluate(as_const(x.view),as_const(residual.view),reductions,ignored);
  okay=okay && status && ignored.accepted;
  LinearConvergenceAuditResult unchanged;unchanged.metric=123;
  status=audit.evaluate(as_const(x.view),as_const(residual.view),reductions,unchanged);
  okay=okay && !status && status.detail==17857 && unchanged.metric==123;
  int all=okay;MPI_Allreduce(MPI_IN_PLACE,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(rank==0)std::cout<<"cold_pressure_audit algorithm="<<unsigned(algorithm)
      <<" baseline_iterations="<<baseline.iterations<<" baseline_metric="<<before.metric
      <<" iterations="<<corrected.iterations<<" rejected="<<corrected.convergence_rejections
      <<" final_metric="<<after.metric<<" passed="<<all<<'\n';
  return all;
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank{},size{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  bool okay=run(LinearAlgorithm::fgmres,rank,size);
  okay=run(LinearAlgorithm::bicgstab,rank,size)&&okay;
  MPI_Finalize();return okay ? 0:1;
}

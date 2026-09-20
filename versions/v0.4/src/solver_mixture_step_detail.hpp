// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_mixture_rows_detail.hpp"
#include <cstdio>

namespace hundun::v04::detail {
struct FrozenScalarProblem {
  const CartesianKernelPlan& kernels;
  const BoundaryPlan& boundary;
  BoundaryStage boundary_stage;
  ConvectionScheme scheme;
  ConstFieldView boundary_velocity;
  const EquationAssemblyContext& context;
  PrimitiveHistory density;
  PrimitiveHistory& scalar;
  FieldView& trial;
  // ESF uses the frozen-density advective form; passive transport uses
  // the conservative density history and divergence form.
  bool frozen_density{true};
};
struct ScalarCorrectionStorage {
  EquationSystemView equation;
  FieldView boundary_variation, rhs, increment, backup;
};
struct ScalarWork {
  double assembly{},closure{},factor{},solve{},refresh{};
  unsigned long long assemblies{},factor_reuses{},operators{},preconditioners{},reductions{};
};
struct ScalarCorrectionRuntime {
  std::vector<ColdPressureRow>& rows;
  ColdPressureDilu& preconditioner;
  HaloEngine& halo;
  SolverWorkspace& workspace;
  ReductionEngine& reductions;
  LinearIdentity identity;
  LinearSolveControl control;
  ScalarWork* work{};
};

inline Status report_scalar_work(const char* module,const ScalarWork& work,MPI_Comm comm) noexcept {
  double local[5]{work.assembly,work.closure,work.factor,work.solve,work.refresh},maximum[5]{};
  unsigned long long count[5]{work.assemblies,work.factor_reuses,work.operators,work.preconditioners,work.reductions},counts[5]{};
  if (MPI_Reduce(local,maximum,5,MPI_DOUBLE,MPI_MAX,0,comm)!=MPI_SUCCESS ||
      MPI_Reduce(count,counts,5,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,comm)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,17866};
  int rank{};
  if (MPI_Comm_rank(comm,&rank)!=MPI_SUCCESS) return {StatusCode::mpi_failure,17866};
  if(rank==0)std::fprintf(stdout,"transport_cost module=%s assembly=%.9g closure=%.9g factor=%.9g solve=%.9g refresh=%.9g assemblies=%llu factor_reuses=%llu operators=%llu preconditioners=%llu reductions=%llu scope=max_rank_stages\n",
      module,maximum[0],maximum[1],maximum[2],maximum[3],maximum[4],counts[0],counts[1],counts[2],counts[3],counts[4]);
  return {};
}

// One correction from the common equation assembly. assemble(certificate)
// reads the current scalar/history and fills storage.equation; refresh(trial)
// closes its boundary/halo after the correction. Caller-owned accepted and
// previous histories stay fixed across repeated corrections. Scratch and the
// DILU factor storage are prepared once by the enclosing compiled scheduler.
template<class Assemble,class Refresh>
LinearSolveResult correct_scalar(const FrozenScalarProblem& p,
    const ScalarCorrectionStorage& storage, ScalarCorrectionRuntime& runtime,
    Assemble&& assemble, Refresh&& refresh,
    const EquationAssemblyCertificate* assembled=nullptr) noexcept {
  constexpr std::uint32_t invalid=17864,numerical=17865;
  LinearSolveResult result;
  const auto cells=p.kernels.cells();
  const auto& q=p.trial;
  const auto ghosts=std::max({q.ghosts.x,q.ghosts.y,q.ghosts.z});
  auto status=valid_cell_view(as_const(q),cells,0,1,1) && q.components==1 &&
      q.revision<UINT64_MAX && ghosts<=UINT8_MAX &&
      p.scalar.trial.base==q.base && p.scalar.trial.field==q.field &&
      p.scalar.trial.revision==q.revision &&
      valid_cell_view(as_const(storage.backup),cells,0,1,std::uint8_t(ghosts)) &&
      valid_cell_view(as_const(storage.rhs),cells,0,1,0) &&
      valid_cell_view(as_const(storage.increment),cells,0,1,1) &&
      runtime.rows.size()==std::size_t(cells.x)*cells.y*cells.z
      ? Status{} : Status{StatusCode::invalid_plan,invalid};
  const std::array<FieldView,7> outputs{storage.equation.diagonal,
      storage.equation.rhs,storage.equation.residual,storage.boundary_variation,
      storage.rhs,storage.increment,storage.backup};
  const std::array<ConstFieldView,7> inputs{as_const(q),p.scalar.accepted,
      p.scalar.previous,p.density.trial,p.density.accepted,p.density.previous,
      p.boundary_velocity};
  for(std::size_t i=1;i<inputs.size();++i)
    if(field_views_overlap(as_const(q),inputs[i])) status={StatusCode::invalid_plan,invalid};
  for(std::size_t i=0;i<outputs.size() && status;++i) {
    if(runtime.workspace.overlaps_storage(outputs[i]) ||
       cell_input_aliases_faces(as_const(outputs[i]),storage.equation,true) ||
       cell_face_views_overlap(as_const(outputs[i]),p.context.mass_flux.x) ||
       cell_face_views_overlap(as_const(outputs[i]),p.context.mass_flux.y) ||
       cell_face_views_overlap(as_const(outputs[i]),p.context.mass_flux.z))
      status={StatusCode::invalid_plan,invalid};
    for(std::size_t j=0;j<i;++j)
      if(field_views_overlap(as_const(outputs[i]),as_const(outputs[j])))
        status={StatusCode::invalid_plan,invalid};
    for(auto input:inputs)
      if(field_views_overlap(input,as_const(outputs[i])))
        status={StatusCode::invalid_plan,invalid};
  }
  for(auto input:inputs)
    if(runtime.workspace.overlaps_storage(input)) status={StatusCode::invalid_plan,invalid};
  status=runtime.reductions.consensus(status);
  if(!status) {
    result.status=status;
    result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();return result;
  }
  const auto original=p.trial;
  const auto original_scalar=p.scalar.trial;
  const auto copy=[&](ConstFieldView input,FieldView output) noexcept {
    for(int z=-q.ghosts.z;z<cells.z+q.ghosts.z;++z)
      for(int y=-q.ghosts.y;y<cells.y+q.ghosts.y;++y)
        for(int x=-q.ghosts.x;x<cells.x+q.ghosts.x;++x)
          output.unchecked({x,y,z},0)=input.unchecked({x,y,z},0);
  };
  copy(as_const(q),storage.backup);
  const auto rollback=[&]() noexcept {
    p.trial=original;p.scalar.trial=original_scalar;
    copy(as_const(storage.backup),p.trial);
  };
  EquationAssemblyCertificate certificate;
  auto clock=MPI_Wtime();
  if(assembled)certificate=*assembled;else {
    status=assemble(certificate);
    if(runtime.work)++runtime.work->assemblies;
  }
  if(runtime.work)runtime.work->assembly+=MPI_Wtime()-clock;
  clock=MPI_Wtime();
  if(status) {
    if(p.frozen_density) status=close_frozen_density_scalar_rows(p.kernels,p.boundary,
        p.boundary_stage,q.field,p.scheme,p.boundary_velocity,p.context,p.density,
        p.scalar.trial,storage.equation,certificate,storage.boundary_variation,
        {runtime.rows.data(),runtime.rows.size()});
    else status=close_mixture_scalar_rows(p.kernels,p.boundary,p.boundary_stage,
        q.field,p.scheme,p.boundary_velocity,p.context,storage.equation,certificate,
        storage.boundary_variation,{runtime.rows.data(),runtime.rows.size()});
  }
  status=runtime.reductions.consensus(status);
  if(!status) {
    rollback();result.status=status;result.termination=LinearTermination::operator_failure;
    result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();return result;
  }
  if(runtime.work)runtime.work->closure+=MPI_Wtime()-clock;
  clock=MPI_Wtime();
  // Row closure already divides the integrated equation by cell volume.
  status=runtime.preconditioner.prepare();
  if(runtime.work) {
    runtime.work->factor+=MPI_Wtime()-clock;
    runtime.work->factor_reuses+=runtime.preconditioner.reused_last_prepare();
  }
  status=runtime.reductions.consensus(status);
  if(!status) {
    rollback();result.status=status;result.termination=LinearTermination::preconditioner_failure;
    result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();return result;
  }
  std::size_t i{};
  for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x,++i) {
    storage.rhs.unchecked({x,y,z},0)=runtime.rows[i].rhs;
    storage.increment.unchecked({x,y,z},0)=0;
  }
  ColdPressureOperator op(runtime.rows,cells,runtime.halo,runtime.identity);
  const LinearSolveInvocation invocation{as_const(storage.rhs),storage.increment,
      runtime.identity,runtime.control};
  clock=MPI_Wtime();
  result=runtime.control.maximum_norm
      ? solve_bicgstab(op,runtime.preconditioner,invocation,runtime.workspace,runtime.reductions)
      : solve_fgmres(op,runtime.preconditioner,invocation,runtime.workspace,runtime.reductions);
  if(runtime.work) {
    runtime.work->solve+=MPI_Wtime()-clock;
    runtime.work->operators+=result.operator_applies;
    runtime.work->preconditioners+=result.preconditioner_applies;
    runtime.work->reductions+=result.reduction_calls;
  }
  status=runtime.reductions.consensus(result.status);
  if(!status) {
    rollback();result.status=status;
    if(result.lowest_failing_rank<0)
      result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();
    return result;
  }
  for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x)
    if(!std::isfinite(q.unchecked({x,y,z},0)+storage.increment.unchecked({x,y,z},0)))
      status={StatusCode::numerical_failure,numerical};
  status=runtime.reductions.consensus(status);
  if(!status) {
    rollback();result.status=status;result.termination=LinearTermination::non_finite;
    result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();return result;
  }
  for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x)
    p.trial.unchecked({x,y,z},0)+=storage.increment.unchecked({x,y,z},0);
  ++p.trial.revision;
  clock=MPI_Wtime();
  status=refresh(p.trial);
  if(runtime.work)runtime.work->refresh+=MPI_Wtime()-clock;
  status=runtime.reductions.consensus(status);
  if(!status) {
    rollback();result.status=status;result.termination=LinearTermination::operator_failure;
    result.lowest_failing_rank=runtime.reductions.lowest_failing_rank();return result;
  }
  p.scalar.trial=as_const(p.trial);
  return result;
}
template<class Assemble,class Refresh>
LinearSolveResult correct_frozen_scalar(const FrozenScalarProblem& p,
    const ScalarCorrectionStorage& storage, ScalarCorrectionRuntime& runtime,
    Assemble&& assemble, Refresh&& refresh) noexcept {
  return correct_scalar(p,storage,runtime,std::forward<Assemble>(assemble),
      std::forward<Refresh>(refresh));
}
} // namespace hundun::v04::detail

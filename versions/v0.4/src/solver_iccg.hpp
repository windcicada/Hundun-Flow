// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_cold.hpp"

namespace hundun::v04::detail {

// REFERENCE cgsol scales each per-volume row and its RHS before IC/CG. Keep the
// caller's original row for its original-equation audit and correction flux.
inline Status scale_iccg_row(const ColdPressureRow& row, double volume,
                             ColdPressureRow& out) noexcept {
  if (!(volume > 0) || !std::isfinite(volume))
    return {StatusCode::invalid_plan, 17866};
  ColdPressureRow scaled = row;
  scaled.diagonal *= volume;
  scaled.rhs *= volume;
  for (auto& a : scaled.neighbour) a *= volume;
  if (!std::isfinite(scaled.diagonal) || !(scaled.diagonal > 0) ||
      !std::isfinite(scaled.rhs) ||
      std::any_of(scaled.neighbour.begin(), scaled.neighbour.end(),
                  [](double a) { return !std::isfinite(a); }))
    return {StatusCode::invalid_plan, 17866};
  out = scaled;
  return {};
}

// Stable global L2 in the original equation units. This retains the same
// absolute/relative stopping scale when IC/CG operates on volume-scaled rows.
inline Status iccg_original_l2(ConstFieldView values, ReductionEngine& reductions,
                               double& norm) noexcept {
  const auto n=values.interior;
  Status local;
  if(!valid_cell_view(values,n,0,1,0))local={StatusCode::invalid_plan,17873};
  auto status=reductions.consensus(local);if(!status)return status;
  double maximum=0, global_maximum=0;
  for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
    const double value=values.unchecked({x,y,z},0);
    if(!std::isfinite(value))local={StatusCode::numerical_failure,17873};
    else maximum=std::max(maximum,std::abs(value));
  }
  status=reductions.checked_max({&maximum,1},{&global_maximum,1},local);
  if(!status)return status;
  if(global_maximum==0) {norm=0;return {};}
  long double sum=0;
  for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
    const double value=values.unchecked({x,y,z},0)/global_maximum;
    sum+=value*value;
  }
  double partial=static_cast<double>(sum),global_sum=0;
  status=reductions.checked_sum({&partial,1},{&global_sum,1});if(!status)return status;
  const double result=global_maximum*std::sqrt(global_sum);
  if(!std::isfinite(result))return {StatusCode::numerical_failure,17873};
  norm=result;return {};
}

class IccgOriginalEquationAudit final : public LinearConvergenceAudit {
 public:
  IccgOriginalEquationAudit(LinearConvergenceAudit& physical, ConstFieldView volumes,
      FieldView residual, double original_l2_limit, PlanFingerprint identity)
      : physical_(physical), volumes_(volumes), residual_(residual),
        limit_(original_l2_limit), identity_(identity) {}
  LinearConvergenceAuditCertificate certificate() const noexcept override {
    return {identity_};
  }
  double last_original_l2() const noexcept {return last_norm_;}
  double original_l2_limit() const noexcept {return limit_;}
  Status evaluate(ConstFieldView solution,ConstFieldView scaled_residual,
      ReductionEngine& reductions,LinearConvergenceAuditResult& out) noexcept override {
    const auto n=volumes_.interior;
    Status local;
    if(!(limit_>0) || !std::isfinite(limit_) || !identity_ ||
        !valid_cell_view(volumes_,n,0,1,0) ||
        !valid_cell_view(scaled_residual,n,0,1,0) ||
        !valid_cell_view(as_const(residual_),n,0,1,0))
      local={StatusCode::invalid_plan,17873};
    auto status=reductions.consensus(local);if(!status)return status;
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      const Int3 c{x,y,z};const double volume=volumes_.unchecked(c,0);
      if(!(volume>0) || !std::isfinite(volume))local={StatusCode::invalid_plan,17873};
      else residual_.unchecked(c,0)=scaled_residual.unchecked(c,0)/volume;
    }
    status=reductions.consensus(local);if(!status)return status;
    double norm=0;
    status=iccg_original_l2(as_const(residual_),reductions,norm);if(!status)return status;
    LinearConvergenceAuditResult physical;
    status=physical_.evaluate(solution,as_const(residual_),reductions,physical);
    if(!status)return status;
    if(!(physical.limit>0) || !std::isfinite(physical.limit) ||
        physical.metric<0 || !std::isfinite(physical.metric) ||
        physical.accepted!=(physical.metric<=physical.limit))
      return {StatusCode::invalid_plan,17873};
    auto result=physical;
    result.metric=std::max(physical.metric/physical.limit,norm/limit_);
    if(!std::isfinite(result.metric))return {StatusCode::numerical_failure,17873};
    result.unscaled_metric=result.metric;
    result.limit=1.;result.accepted=physical.accepted && norm<=limit_;
    last_norm_=norm;out=result;return {};
  }
 private:
  LinearConvergenceAudit& physical_;
  ConstFieldView volumes_;
  FieldView residual_;
  double limit_,last_norm_{};
  PlanFingerprint identity_;
};

struct IccgMatrixAdmissionReport {
  std::uint32_t grounding_exchanges{};
  bool spd{};
};

// Sufficient global SPD proof for the volume-scaled pressure M matrix:
// exact opposite-face symmetry, weak row dominance, and one strictly
// dominant row in every connected component. Local union/find collapses
// each rank before propagating grounded flags over the actual halo graph.
// Exact symmetry keeps admission distinct from changing the matrix.
// exchange owns the supplied scratch halo; all physical exterior ghosts
// start at zero. The final collective result is uniform across ranks.
// Storage allocations are made by callers before entering this collective
// routine so allocation failures can participate in ordinary consensus.
template<class Exchange>
Status admit_iccg_matrix(const std::vector<ColdPressureRow>& rows,
    MeshPatch patch, Int3 global_cells, const std::array<bool,3>& periodic,
    FieldView scratch, Span<std::size_t> parent, Span<std::uint8_t> grounded,
    Exchange&& exchange, ReductionEngine& reductions,
    IccgMatrixAdmissionReport& report) noexcept {
  report = {};
  const auto n=patch.cells;
  const auto count=std::size_t(n.x)*n.y*n.z;
  Status local;
  if(n.x<=0 || n.y<=0 || n.z<=0 || rows.size()!=count ||
      parent.size!=count || grounded.size!=count || !parent.data || !grounded.data ||
      !valid_cell_view(as_const(scratch),n,0,1,1))
    local={StatusCode::invalid_plan,17869};
  const int starts[]{patch.begin.x,patch.begin.y,patch.begin.z};
  const int extents[]{global_cells.x,global_cells.y,global_cells.z};
  const int widths[]{n.x,n.y,n.z};
  for(unsigned a=0;a<3;++a)
    if(starts[a]<0 || extents[a]<=0 || std::int64_t(starts[a])+widths[a]>extents[a])
      local={StatusCode::invalid_plan,17869};
  auto status=reductions.consensus(local);if(!status)return status;
  const auto flat=[&](Int3 c){return std::size_t(c.x)+n.x*(c.y+std::size_t(n.y)*c.z);};
  const auto root=[&](std::size_t i) {
    while(parent.data[i]!=i) {parent.data[i]=parent.data[parent.data[i]];i=parent.data[i];}
    return i;
  };
  for(std::size_t i=0;i<count;++i) {parent.data[i]=i;grounded.data[i]=0;}
  for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
    const Int3 c{x,y,z};const auto i=flat(c);const auto& row=rows[i];
    long double sum=0;
    if(!(row.diagonal>0) || !std::isfinite(row.diagonal))local={StatusCode::invalid_plan,17869};
    for(unsigned f=0;f<6;++f) {
      const double a=row.neighbour[f];
      if(!(a>=0) || !std::isfinite(a))local={StatusCode::invalid_plan,17869};
      sum+=a;
      Int3 nb=c;auto& coord=f/2==0 ? nb.x : f/2==1 ? nb.y : nb.z;
      coord+=f%2 ? 1:-1;
      const int start=f/2==0 ? patch.begin.x : f/2==1 ? patch.begin.y : patch.begin.z;
      const int extent=f/2==0 ? global_cells.x : f/2==1 ? global_cells.y : global_cells.z;
      if(!periodic[f/2] && (coord+start<0 || coord+start>=extent) && a!=0)
        local={StatusCode::invalid_plan,17869};
      if(f%2 && a>0 && nb.x<n.x && nb.y<n.y && nb.z<n.z)
        parent.data[root(flat(nb))]=root(i);
    }
    if(row.diagonal<sum)local={StatusCode::invalid_plan,17870};
    grounded.data[i]=row.diagonal>sum;
  }
  status=reductions.consensus(local);if(!status)return status;
  const auto clear_ghosts=[&]() {
    for(int z=-1;z<=n.z;++z)for(int y=-1;y<=n.y;++y)for(int x=-1;x<=n.x;++x)
      if(int(x<0 || x>=n.x)+int(y<0 || y>=n.y)+int(z<0 || z>=n.z)==1)
        scratch.unchecked({x,y,z},0)=0;
  };
  for(unsigned axis=0;axis<3;++axis) {
    clear_ghosts();
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x)
      scratch.unchecked({x,y,z},0)=rows[flat({x,y,z})].neighbour[2*axis];
    status=exchange(scratch);status=reductions.consensus(status);if(!status)return status;
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      Int3 nb{x,y,z};(axis==0 ? nb.x : axis==1 ? nb.y : nb.z)++;
      if(rows[flat({x,y,z})].neighbour[2*axis+1]!=scratch.unchecked(nb,0))
        local={StatusCode::invalid_plan,17871};
    }
    status=reductions.consensus(local);if(!status)return status;
  }
  for(std::size_t i=0;i<count;++i)if(grounded.data[i])grounded.data[root(i)]=1;
  const auto ungrounded=[&]() {
    for(std::size_t i=0;i<count;++i)if(!grounded.data[root(i)])return 1.;
    return 0.;
  };
  double missing=ungrounded(), global_missing=0;
  status=reductions.checked_max({&missing,1},{&global_missing,1});if(!status)return status;
  while(global_missing) {
    clear_ghosts();
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x)
      scratch.unchecked({x,y,z},0)=grounded.data[root(flat({x,y,z}))];
    status=exchange(scratch);status=reductions.consensus(status);if(!status)return status;
    ++report.grounding_exchanges;
    double changed=0;
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      const Int3 c{x,y,z};const auto i=flat(c),r=root(i);
      if(grounded.data[r])continue;
      for(unsigned f=0;f<6;++f)if(rows[i].neighbour[f]>0) {
        Int3 nb=c;(f/2==0 ? nb.x : f/2==1 ? nb.y : nb.z)+=f%2 ? 1:-1;
        if(scratch.unchecked(nb,0)==1) {grounded.data[r]=1;changed=1;break;}
      }
    }
    double signals[]{ungrounded(),changed}, global[2]{};
    status=reductions.checked_max({signals,2},{global,2});if(!status)return status;
    global_missing=global[0];
    if(global_missing && !global[1])return {StatusCode::invalid_plan,17872};
    if(report.grounding_exchanges==std::numeric_limits<std::uint32_t>::max())
      return {StatusCode::invalid_plan,17872};
  }
  report.spd=true;
  return {};
}

// Block-local IC(0), matching cgsol's x/y/z triangular ordering. MPI and
// periodic wraps remain in the operator; their ghost factors are zero in
// the original routine. The backward sweep explicitly uses L transpose,
// so positive pivots certify an SPD preconditioner independently of A.
// The operator must separately pass its global SPD admission before PCG.
class ColdPressureIccg final : public LinearPreconditioner {
 public:
  ColdPressureIccg(const std::vector<ColdPressureRow>& rows, Int3 cells,
                   LinearIdentity identity)
      : rows_(rows), cells_(cells), identity_(identity), inverse_(rows.size()) {}
  void reset_identity(LinearIdentity identity) noexcept {identity_=identity;prepared_=false;}
  Status prepare() noexcept {
    prepared_ = false;
    if (cells_.x <= 0 || cells_.y <= 0 || cells_.z <= 0 ||
        rows_.size() != std::size_t(cells_.x)*cells_.y*cells_.z)
      return {StatusCode::invalid_plan, 17867};
    const std::size_t strides[]{1, std::size_t(cells_.x),
                                std::size_t(cells_.x)*cells_.y};
    std::size_t i = 0;
    for (int z=0; z<cells_.z; ++z) for (int y=0; y<cells_.y; ++y)
      for (int x=0; x<cells_.x; ++x, ++i) {
        const int at[]{x,y,z};
        double pivot = rows_[i].diagonal;
        for (unsigned a=0; a<3; ++a) {
          const double lower = rows_[i].neighbour[2*a];
          if (!std::isfinite(lower)) return {StatusCode::invalid_plan,17867};
          if (at[a]>0) pivot -= lower*lower*inverse_[i-strides[a]];
        }
        if (!(pivot>0) || !std::isfinite(pivot) || !std::isfinite(1/pivot))
          return {StatusCode::numerical_failure,17868};
        inverse_[i] = 1/pivot;
      }
    prepared_ = true;
    return {};
  }
  Span<const double> inverse_pivots() const noexcept {
    return prepared_ ? Span<const double>{inverse_.data(),inverse_.size()}
                     : Span<const double>{};
  }
  LinearPreconditionerCertificate certificate() const noexcept override {
    return {prepared_ ? identity_ : LinearIdentity{}, UINT64_C(0x434f4c4449434331),
            LinearPreconditionerClass::fixed_spd};
  }
  Status apply(ConstFieldView input, FieldView output, std::uint32_t) noexcept override {
    if (!prepared_ || !valid_cell_view(input,cells_,0,1,0) ||
        !valid_cell_view(as_const(output),cells_,0,1,0))
      return {StatusCode::invalid_plan,17867};
    std::size_t i=0;
    for (int z=0; z<cells_.z; ++z) for (int y=0; y<cells_.y; ++y)
      for (int x=0; x<cells_.x; ++x,++i) {
        double value=input.unchecked({x,y,z},0);
        if (x>0) value+=rows_[i].neighbour[0]*output.unchecked({x-1,y,z},0);
        if (y>0) value+=rows_[i].neighbour[2]*output.unchecked({x,y-1,z},0);
        if (z>0) value+=rows_[i].neighbour[4]*output.unchecked({x,y,z-1},0);
        output.unchecked({x,y,z},0)=value*inverse_[i];
      }
    const std::size_t sy=cells_.x, sz=sy*cells_.y;
    i=rows_.size();
    for (int z=cells_.z-1; z>=0; --z) for (int y=cells_.y-1; y>=0; --y)
      for (int x=cells_.x-1; x>=0; --x) {
        --i; double value=0;
        if (x+1<cells_.x) value+=rows_[i+1].neighbour[0]*output.unchecked({x+1,y,z},0);
        if (y+1<cells_.y) value+=rows_[i+sy].neighbour[2]*output.unchecked({x,y+1,z},0);
        if (z+1<cells_.z) value+=rows_[i+sz].neighbour[4]*output.unchecked({x,y,z+1},0);
        output.unchecked({x,y,z},0)+=value*inverse_[i];
      }
    return {};
  }
 private:
  const std::vector<ColdPressureRow>& rows_;
  Int3 cells_;
  LinearIdentity identity_;
  std::vector<double> inverse_;
  bool prepared_{};
};
}  // namespace hundun::v04::detail

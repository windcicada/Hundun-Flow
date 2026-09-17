// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_cold.hpp"

namespace hundun::v04::detail {

// COAST cgsol scales each per-volume row and its RHS before IC/CG. Keep the
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

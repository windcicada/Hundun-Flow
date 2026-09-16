// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core_tcr_dynamic_history_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "hundun/v04_ibm.hpp"
#include <climits>

namespace hundun::v04::detail {
// Cold-owned scalar statistics exchange. Source and target identities use
// global Cartesian cells, including corners and odd MPI partitions.
class DynamicTcrPlan {
public:
  Status configure(MPI_Comm comm, const CartesianGeometryPlan &geometry,
      MeshPatch patch, std::array<bool,3> periodic, std::size_t species,
      std::array<std::size_t,3> groups, std::uint64_t maximum_bytes);
  Status finish(DynamicTcrHistory &, Span<const FieldView> fields,
      ConstFieldView auxiliary, ConstFieldView density,
      Span<const std::uint8_t> activity, std::uint64_t step) noexcept;
  std::uint64_t owned_bytes() const noexcept { return owned_bytes_; }
  RemoteDonorExchangeStats stats() const noexcept { return halo_.stats(); }
private:
  Status agree(Status) const noexcept;
  Status exchange() noexcept;
  bool canonical(Int3 local, Int3 &global) const noexcept;
  Int3 neighbor(Int3 cell, unsigned axis, int direction) const noexcept;
  double coordinate(Int3 cell, unsigned axis) const noexcept;
  template<class F> void filter(Int3 c, F consume) const noexcept {
    // A global low face uses a forward pair. Partition faces use the same
    // backward test filter as interior cells.
    int offset[3]{-1,-1,-1};
    const int local[]{c.x,c.y,c.z}, begin[]{patch_.begin.x,patch_.begin.y,patch_.begin.z};
    for(unsigned d=0;d<3;++d)
      if(!periodic_[d] && local[d]+begin[d]==0) offset[d]=1;
    for(int k=0;k<2;++k)for(int j=0;j<2;++j)for(int i=0;i<2;++i) {
      const Int3 p{c.x+i*offset[0],c.y+j*offset[1],c.z+k*offset[2]};
      Int3 g{};
      if(canonical(p,g) && view_.unchecked(p,0)>0) consume(p);
    }
  }
  static Status invalid() noexcept { return {StatusCode::invalid_plan,10237}; }
  MPI_Comm comm_{MPI_COMM_NULL};
  int rank_{};
  const CartesianGeometryPlan *geometry_{};
  MeshPatch patch_{};
  std::array<bool,3> periodic_{};
  std::array<std::size_t,3> groups_{};
  std::size_t ns_{}, count_{};
  std::uint64_t owned_bytes_{};
  std::vector<double> storage_;
  FieldView view_{};
  RemoteDonorExchangePlan halo_;
};
} // namespace hundun::v04::detail

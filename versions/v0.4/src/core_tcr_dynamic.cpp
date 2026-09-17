// SPDX-License-Identifier: Apache-2.0
#include "core_tcr_dynamic_detail.hpp"
#include <limits>

namespace hundun::v04::detail {
Status DynamicTcrPlan::agree(Status local) const noexcept {
  int failed=local ? INT_MAX : rank_, first{};
  if(MPI_Allreduce(&failed,&first,1,MPI_INT,MPI_MIN,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  if(first==INT_MAX)return {};
  unsigned words[]{unsigned(local.code),local.detail};
  if(MPI_Bcast(words,2,MPI_UNSIGNED,first,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  return {static_cast<StatusCode>(words[0]),words[1]};
}
bool DynamicTcrPlan::canonical(Int3 local,Int3 &global) const noexcept {
  const auto n=geometry_->global_cells();
  int p[]{local.x+patch_.begin.x,local.y+patch_.begin.y,local.z+patch_.begin.z};
  const int size[]{n.x,n.y,n.z};
  for(unsigned d=0;d<3;++d) {
    if(periodic_[d])p[d]=(p[d]%size[d]+size[d])%size[d];
    else if(p[d]<0 || p[d]>=size[d])return false;
  }
  global={p[0],p[1],p[2]};return true;
}
Int3 DynamicTcrPlan::neighbor(Int3 p,unsigned d,int direction) const noexcept {
  if(d==0)p.x+=direction;else if(d==1)p.y+=direction;else p.z+=direction;
  return p;
}
double DynamicTcrPlan::coordinate(Int3 p,unsigned d) const noexcept {
  const int g[]{p.x+patch_.begin.x,p.y+patch_.begin.y,p.z+patch_.begin.z};
  const auto &axis=geometry_->axis(static_cast<CartesianAxis>(d));
  const auto centres=axis.centres(),faces=axis.faces();
  const auto n=int(centres.size);
  const int index=(g[d]%n+n)%n;
  const int cycle=(g[d]-index)/n;
  return centres.data[index]+cycle*(faces.data[faces.size-1]-faces.data[0]);
}
Status DynamicTcrPlan::exchange() noexcept {
  auto s=agree(halo_.preflight_exchange(10,{&view_,1}));
  if(s)s=halo_.exchange(10,{&view_,1});
  return s;
}
Status DynamicTcrPlan::configure(MPI_Comm comm,const CartesianGeometryPlan &geometry,
    MeshPatch patch,std::array<bool,3> periodic,std::size_t species,
    std::array<std::size_t,3> groups,std::uint64_t maximum_bytes) {
  if(comm==MPI_COMM_NULL)return invalid();
  comm_=comm;
  if(MPI_Comm_rank(comm,&rank_)!=MPI_SUCCESS)return {StatusCode::mpi_failure,10237};
  geometry_=&geometry;patch_=patch;periodic_=periodic;ns_=species;groups_=groups;
  auto s=agree(ns_<2 || ns_>120 || !geometry.fingerprint() ? invalid() : Status{});
  if(!s)return s;
  const std::uint64_t sy=std::uint64_t(patch.cells.x)+4,
      sz=sy*(std::uint64_t(patch.cells.y)+4),stride=sz*(std::uint64_t(patch.cells.z)+4);
  const auto components=std::max<std::size_t>(20,1+2*ns_);
  s=agree(stride>SIZE_MAX/(8*components) || stride>maximum_bytes/(8*components)
      ? Status{StatusCode::allocation_failure,10237} : Status{});
  if(!s)return s;
  std::vector<GlobalCellId> ids;std::vector<Int3> targets;
  try {
    storage_.assign(stride*components,0.);
    for(int z=-2;z<patch.cells.z+2;++z)for(int y=-2;y<patch.cells.y+2;++y)
      for(int x=-2;x<patch.cells.x+2;++x) {
        if(x>=0 && x<patch.cells.x && y>=0 && y<patch.cells.y && z>=0 && z<patch.cells.z)continue;
        const Int3 p{x,y,z};Int3 g{};
        if(!canonical(p,g))continue;
        const auto n=geometry.global_cells();
        ids.push_back(std::uint64_t(g.x)+std::uint64_t(n.x)*(std::uint64_t(g.y)+std::uint64_t(n.y)*g.z));
        targets.push_back(p);
      }
  } catch(...) { s={StatusCode::allocation_failure,10237}; }
  s=agree(s);if(!s)return s;
  view_.base=storage_.data()+2+2*sy+2*sz;view_.interior=patch.cells;view_.ghosts={2,2,2};
  view_.components=static_cast<std::uint8_t>(components);view_.stride_y=sy;view_.stride_z=sz;
  view_.component_stride=stride;view_.field=1;view_.revision=1;
  view_.storage_identity=reinterpret_cast<StorageIdentity>(storage_.data());
  view_.revision_domain=geometry.fingerprint();
  const RemoteDonorTargets target{geometry.fingerprint(),2,periodic,{ids.data(),ids.size()},{targets.data(),targets.size()}};
  const RemoteDonorFieldSpec field{1,static_cast<std::uint8_t>(components)};
  s=RemoteDonorExchangePlan::analyze_cells(comm,geometry.global_cells(),patch,target,{&field,1},10,halo_);
  if(!s)return s;
  const auto stats=halo_.stats();
  owned_bytes_=sizeof(*this)+8*storage_.capacity()+256*(stats.received_cells+stats.supplied_cells)+2*stats.bytes_per_exchange;
  s=agree(owned_bytes_>maximum_bytes ? Status{StatusCode::allocation_failure,10237} : Status{});
  if(s)s=halo_.bind(comm);
  count_=std::size_t(patch.cells.x)*patch.cells.y*patch.cells.z;
  return s;
}
Status DynamicTcrPlan::stage_rates(Dyn711History &history,
    const TurbulencePlan &turbulence, const Dyn711StatisticsInput &input) noexcept {
  using namespace tcr::detail;
  if(comm_==MPI_COMM_NULL)return invalid();
  const auto cells=patch_.cells;
  auto status=agree(history.cells()!=count_ || history.species()!=ns_ ||
      !std::isfinite(input.dt) || input.dt<=0 ||
      !std::isfinite(input.weak_rate) || input.weak_rate<=0 ||
      input.molecular_weights.size!=ns_ || !input.molecular_weights.data ||
      !valid_cell_view(input.physical_mass_fractions,cells,0,ns_,0) ||
      !valid_cell_view(input.pdf_rates,cells,0,ns_,0) ||
      !valid_cell_view(input.psr_rates,cells,0,ns_,0) ||
      !valid_cell_view(input.eta,cells,0,1,0) ||
      !valid_cell_view(input.density,cells,0,1,0) ||
      !valid_cell_view(input.molecular_viscosity,cells,0,1,0) ||
      !valid_cell_view(input.velocity_gradient,cells,0,9,0) ||
      (input.activity.size && (input.activity.size!=count_ || !input.activity.data))
      ? invalid() : Status{});
  if(!status)return status;
  // Before any conditional model work, ensure each rank advances the same
  // accepted model clock and interval. Fixed-size reductions precede ns data.
  std::uint64_t lower[2]{history.statistics_calls(),ns_},upper[2]{lower[0],lower[1]};
  if(MPI_Allreduce(MPI_IN_PLACE,lower,2,MPI_UINT64_T,MPI_MIN,comm_)!=MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE,upper,2,MPI_UINT64_T,MPI_MAX,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  double controls[4]{input.dt,input.weak_rate,-input.dt,-input.weak_rate};
  if(MPI_Allreduce(MPI_IN_PLACE,controls,4,MPI_DOUBLE,MPI_MAX,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  if(lower[0]!=upper[0] || lower[1]!=upper[1] ||
      controls[0]!=-controls[2] || controls[1]!=-controls[3])
    return invalid();
  for(std::size_t q=0;q<ns_;++q)
    if(!std::isfinite(input.molecular_weights.data[q]) || input.molecular_weights.data[q]<=0)
      status=invalid();
  status=agree(status);if(!status)return status;
  const auto active=[&](std::size_t i){return !input.activity.size || input.activity.data[i]!=0;};
  // Floors and the one-second fallback cap are the original chemical-clock
  // policy. The physical local SGS clock has its independently approved units.
  std::array<double,240> maxima;
  maxima.fill(1e-12);
  const auto finite_time=[](long double value) {
    return static_cast<double>(std::min(value,
        static_cast<long double>(std::numeric_limits<double>::max())));
  };
  std::size_t i{};
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(input.activity.size && input.activity.data[i]>1)status=invalid();
    if(!active(i))continue;
    const Int3 cell{x,y,z};
    const double eta=input.eta.unchecked(cell,0);
    if(!std::isfinite(eta) || eta<0 || eta>1)status=invalid();
    VelocityGradient gradient;
    for(unsigned q=0;q<9;++q)gradient.value[q]=input.velocity_gradient.unchecked(cell,q);
    SgsState sgs;
    const double rho=input.density.unchecked(cell,0),mu=input.molecular_viscosity.unchecked(cell,0);
    const auto evaluated=turbulence.evaluate_sgs_cell(cell,gradient,rho,mu,sgs);
    if(!evaluated) {status=evaluated;continue;}
    const auto times=dyn711_flow_times(sgs.kinetic_energy_m2_s2,sgs.dissipation_w_m3,rho,mu);
    if(!times.available) {status=invalid();continue;}
    view_.unchecked(cell,0)=times.flow;
    double sum{};
    for(std::size_t q=0;q<ns_;++q) {
      const double yq=input.physical_mass_fractions.unchecked(cell,q);
      const double rate=input.psr_rates.unchecked(cell,q),pdf=input.pdf_rates.unchecked(cell,q);
      if(!std::isfinite(yq) || yq<0 || yq>1 || !std::isfinite(rate) || !std::isfinite(pdf)) {
        status=invalid();continue;
      }
      sum+=yq;
      const double n=yq/input.molecular_weights.data[q];
      if(!std::isfinite(n)) {status=invalid();continue;}
      maxima[ns_+q]=std::max(maxima[ns_+q],n);
      if(std::abs(rate)>=input.weak_rate && n>=input.weak_rate)
        maxima[q]=std::max(maxima[q],finite_time(static_cast<long double>(n)/std::abs(rate)));
    }
    if(std::abs(sum-1)>2e-12)status=invalid();
  }
  status=agree(status);if(!status)return status;
  if(MPI_Allreduce(MPI_IN_PLACE,maxima.data(),static_cast<int>(2*ns_),MPI_DOUBLE,MPI_MAX,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(!active(i)) {
      const auto staged=history.stage_inactive(i);if(!staged)status=staged;
      continue;
    }
    const Int3 cell{x,y,z};
    for(std::size_t q=0;q<ns_;++q) {
      const double rate=input.psr_rates.unchecked(cell,q);
      const double n=input.physical_mass_fractions.unchecked(cell,q)/input.molecular_weights.data[q];
      const double chemical=std::abs(rate)<input.weak_rate || n<input.weak_rate
          ? std::min(1.,maxima[q]) : finite_time(rate>0
              ? .5L*maxima[ns_+q]/rate : static_cast<long double>(n)/std::abs(rate));
      const auto staged=history.stage_rate(i,q,input.dt,input.pdf_rates.unchecked(cell,q),rate,
          input.eta.unchecked(cell,0),chemical,view_.unchecked(cell,0),input.weak_rate);
      if(!staged)status=staged;
    }
  }
  return agree(status);
}
Status DynamicTcrPlan::finish(DynamicTcrHistory &history,Span<const FieldView> fields,
    ConstFieldView auxiliary,ConstFieldView density,Span<const std::uint8_t> activity,
    std::uint64_t step) noexcept {
  using namespace tcr::detail;
  const auto cells=patch_.cells;
  auto s=agree(history.species()!=ns_ || !fields.size || !fields.data ||
      !valid_cell_view(auxiliary,cells,0,ns_+1,0) || !valid_cell_view(density,cells,0,1,0) ||
      (activity.size && (activity.size!=count_ || !activity.data)) ? invalid() : Status{});
  if(!s)return s;
  for(std::size_t f=0;f<fields.size;++f)
    if(!valid_cell_view(as_const(fields.data[f]),cells,0,ns_+1,0))s=invalid();
  s=agree(s);if(!s)return s;
  const auto active=[&](std::size_t i){return !activity.size || activity.data[i]!=0;};
  // Complete weak species rates from the same species' eligible neighbors,
  // with the declared all-species global mean as the COAST fallback.
  double totals[2]{};std::size_t i{};
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    const Int3 c{x,y,z};const auto *p=history.candidate(i);
    view_.unchecked(c,0)=active(i) ? 1 : 0;
    for(std::size_t q=0;q<ns_;++q) {
      view_.unchecked(c,1+2*q)=p[5*q+3];view_.unchecked(c,2+2*q)=p[5*q+4];
      if(active(i) && (p[5*q+4]==1 || p[5*q+4]==2)) {totals[0]+=p[5*q+3];totals[1]+=1;}
    }
  }
  if(MPI_Allreduce(MPI_IN_PLACE,totals,2,MPI_DOUBLE,MPI_SUM,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  s=exchange();if(!s)return s;
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(!active(i))continue;
    auto *p=history.candidate(i);
    for(std::size_t q=0;q<ns_;++q) {
      if(p[5*q+4]!=3)continue;
      double sum{},count{};
      for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx) {
        const Int3 c{x+dx,y+dy,z+dz};Int3 g{};
        if(!canonical(c,g) || view_.unchecked(c,0)==0)continue;
        const auto state=view_.unchecked(c,2+2*q);
        if(state==1 || state==2){sum+=view_.unchecked(c,1+2*q);count+=1;}
      }
      const double value=count ? sum/count : totals[1] ? totals[0]/totals[1] : 1.;
      p[5*q+2]=p[5*q+3]=std::clamp(value,1e-4,1.);
    }
  }
  if(step%4!=0)return history.seal();
  // rho, V, physical means (3), field0 coordinates (3), gradients (9),
  // and two products per group share one prepared statistics arena.
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    const Int3 c{x,y,z};const double rho=density.unchecked(c,0);
    if(active(i) && (!std::isfinite(rho) || rho<=0))s=invalid();
    view_.unchecked(c,0)=active(i) ? rho : 0.;
    view_.unchecked(c,1)=geometry_->x().widths().data[x+patch_.begin.x]*
        geometry_->y().widths().data[y+patch_.begin.y]*geometry_->z().widths().data[z+patch_.begin.z];
    for(unsigned group=0;group<3;++group) {
      double mean{};
      if(groups_[group]<ns_)for(std::size_t f=0;f<fields.size;++f)
        mean+=fields.data[f].unchecked(c,groups_[group])/fields.size;
      view_.unchecked(c,2+group)=mean;
      view_.unchecked(c,5+group)=groups_[group]<ns_ ? auxiliary.unchecked(c,groups_[group]) : 0.;
    }
  }
  s=agree(s);if(s)s=exchange();if(!s)return s;
  s=spatial_products(false);if(!s)return s;
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(!active(i))continue;
    auto *p=history.candidate(i);
    for(unsigned group=0;group<3;++group) {
      if(groups_[group]>=ns_)continue;
      double m2{},lm{};
      filter({x,y,z},[&](Int3 q){m2+=view_.unchecked(q,17+group);lm+=view_.unchecked(q,5+group);});
      p[5*ns_+group]=dynamic_cd_from_products(m2,lm);
    }
  }
  return history.seal();
}

Status DynamicTcrPlan::spatial_products(bool dyn711) noexcept {
  using namespace tcr::detail;
  const auto cells=patch_.cells;
  Status s;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
    const Int3 c{x,y,z};
    for(unsigned group=0;group<3;++group)for(unsigned d=0;d<3;++d) {
      auto a=neighbor(c,d,-1),b=neighbor(c,d,1);Int3 g{};
      if(!canonical(a,g) || view_.unchecked(a,0)==0)a=c;
      if(!canonical(b,g) || view_.unchecked(b,0)==0)b=c;
      const double distance=coordinate(b,d)-coordinate(a,d);
      view_.unchecked(c,8+3*group+d)=distance>0 ?
          (view_.unchecked(b,5+group)-view_.unchecked(a,5+group))/distance : 0.;
    }
  }
  s=exchange();if(!s)return s;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
    const Int3 c{x,y,z};
    for(unsigned group=0;group<3;++group) {
      std::array<DynamicFilterDonor,8> donors;
      unsigned count{};
      filter(c,[&](Int3 p) {
        auto &donor=donors[count++];
        donor.density=view_.unchecked(p,0);
        donor.volume=view_.unchecked(p,1);
        donor.scalar=view_.unchecked(p,2+group);
        for(unsigned d=0;d<3;++d)donor.gradient[d]=view_.unchecked(p,8+3*group+d);
      },dyn711);
      DynamicFilterProducts products;
      if(count) {
        DynamicFilterMoments moments;
        if(!dynamic_filter_moments(donors.data(),count,moments))s=invalid();
        else {
          products=dyn711 ? dyn711_filter_products(moments) : dynamic_filter_products(moments);
          if(!products.available)s=invalid();
        }
      }
      // Use the vacated field0/physical coordinate slots after all gradient
      // reads. Separate products storage avoids aliasing neighboring filters.
      view_.unchecked(c,17+group)=products.available ? products.m_squared : 0.;
      view_.unchecked(c,5+group)=products.available ? products.l_times_m : 0.;
    }
  }
  s=agree(s);if(s)s=exchange();if(!s)return s;
  return s;
}

Status DynamicTcrPlan::finish(Dyn711History &history,
    const std::array<ConstFieldView,3> &means,
    const std::array<ConstFieldView,3> &gradient_coordinates,
    ConstFieldView density,Span<const std::uint8_t> activity) noexcept {
  using namespace tcr::detail;
  const auto cells=patch_.cells;
  auto s=agree(history.species()!=ns_ || history.cells()!=count_ ||
      !valid_cell_view(density,cells,0,1,0) ||
      (activity.size && (activity.size!=count_ || !activity.data)) ? invalid() : Status{});
  if(!s)return s;
  for(unsigned g=0;g<3;++g)
    if(!valid_cell_view(means[g],cells,0,1,0) ||
        !valid_cell_view(gradient_coordinates[g],cells,0,1,0))s=invalid();
  s=agree(s);if(!s)return s;
  // The model clock controls collective filter exchanges. Check it across
  // partitions before taking that branch, including imported/rebuilt history.
  const std::uint64_t calls=history.statistics_calls();
  std::uint64_t minimum{},maximum{};
  if(MPI_Allreduce(&calls,&minimum,1,MPI_UINT64_T,MPI_MIN,comm_)!=MPI_SUCCESS ||
      MPI_Allreduce(&calls,&maximum,1,MPI_UINT64_T,MPI_MAX,comm_)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,10237};
  if(minimum!=maximum)return invalid();
  if(!dyn711_tick(history.clock()).update_cphi)return agree(history.seal());
  const auto active=[&](std::size_t i){return !activity.size || activity.data[i]!=0;};
  std::size_t i{};
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    const Int3 c{x,y,z};const double rho=density.unchecked(c,0);
    if(active(i) && (!std::isfinite(rho) || rho<=0))s=invalid();
    view_.unchecked(c,0)=active(i) ? rho : 0.;
    view_.unchecked(c,1)=geometry_->x().widths().data[x+patch_.begin.x]*
        geometry_->y().widths().data[y+patch_.begin.y]*geometry_->z().widths().data[z+patch_.begin.z];
    for(unsigned g=0;g<3;++g) {
      view_.unchecked(c,2+g)=means[g].unchecked(c,0);
      view_.unchecked(c,5+g)=gradient_coordinates[g].unchecked(c,0);
    }
  }
  s=agree(s);if(s)s=exchange();if(!s)return s;
  s=spatial_products(true);if(!s)return s;
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(!active(i))continue;
    const Int3 c{x,y,z};std::array<double,3> ratios;
    for(unsigned g=0;g<3;++g) {
      double m2{},lm{};
      filter(c,[&](Int3 p){m2+=view_.unchecked(p,17+g);lm+=view_.unchecked(p,5+g);},true);
      ratios[g]=dyn711_filter_ratio(m2,lm);
    }
    // Coordinates have been consumed; a separate component avoids overwriting
    // products still read by adjacent second filters.
    const double value=dyn711_select_cphi(ratios);
    if(!std::isfinite(value))s=invalid();
    view_.unchecked(c,2)=value;
  }
  s=agree(s);if(s)s=exchange();if(!s)return s;
  i=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
    if(!active(i))continue;
    const Int3 c{x,y,z};std::array<double,7> neighbors;
    // Missing boundary/solid donors have zero smoothing weight. A value at
    // the endpoint 1 is outside the source's strict (1.5,12) inclusion band.
    neighbors.fill(1.);unsigned count{};
    filter(c,[&](Int3 p){if(p.x!=x || p.y!=y || p.z!=z)neighbors[count++]=view_.unchecked(p,2);},true);
    const auto staged=history.stage_cphi(i,dyn711_smooth_cphi(view_.unchecked(c,2),neighbors));
    if(!staged)s=staged;
  }
  s=agree(s);if(s)s=agree(history.seal());return s;
}

} // namespace hundun::v04::detail

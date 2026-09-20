// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include "solver_equation_detail.hpp"
#include "solver_scalar_boundary_detail.hpp"
#include <limits>

namespace hundun::v04::detail {

// Attempt-local physical ensemble ledger. Sources freeze once; each pressure
// solve replaces its correction. Only the enclosing accepted transaction
// publishes the report. The dependent species follows mixture mass closure.
class CompositionBalanceLedger {
 public:
  enum Term : unsigned { accepted, current, temporal, transport, pressure,
                         noise, mixing, chemistry, term_count };
  using Row = std::array<long double, term_count>;
  Status initialize(const portable::GasIdentity& identity,
      Span<const std::size_t> mapping, std::size_t dependent,
      std::size_t fields, RevisionToken generation, double dt,
      bool after_parcels) noexcept {
    const auto ns = identity.species_names.size();
    const auto ne = identity.element_names.size();
    if (ns < 2 || ns >= UINT8_MAX || !ne || ne >= UINT8_MAX ||
        mapping.size + 1 != ns || !mapping.data || dependent >= ns ||
        !fields || !generation || !(dt > 0) || !std::isfinite(dt) ||
        identity.molecular_weights_kg_per_kmol.size() != ns ||
        identity.element_counts.size() != ns * ne) return invalid();
    std::array<bool, UINT8_MAX> seen{};
    seen[dependent] = true;
    for (std::size_t i=0;i<mapping.size;++i) {
      const auto s=mapping.data[i];
      if (s >= ns || seen[s]) return invalid();
      seen[s] = true;
    }
    for (const auto mw : identity.molecular_weights_kg_per_kmol)
      if (!(mw > 0) || !std::isfinite(mw)) return invalid();
    try {
      rows_.assign(ns + 1, {});
      storage_bounds_.assign(ns + 1, 0);
      mapping_.assign(mapping.data, mapping.data + mapping.size);
      transport_counts_.assign(ns, 0);
      pressure_counts_.assign(ns, 0);
    } catch (...) { return {StatusCode::allocation_failure, 17861}; }
    identity_ = &identity; dependent_ = dependent; fields_ = fields;
    generation_ = generation; dt_ = dt; after_parcels_ = after_parcels;
    return {};
  }
  bool active() const noexcept { return identity_ != nullptr; }
  void add(std::size_t independent, Term term, long double value) noexcept {
    rows_[mapping_[independent]][term] += value;
  }
  void add_total(Term term, long double value) noexcept { rows_.back()[term] += value; }
  // Each local storage product carries the resolution of its FP64 inputs.
  // eps*abs(x)+denorm_min encloses one representable spacing, including zero.
  // The ensemble mean has the mean input uncertainty; physical Y >= 0 makes
  // that envelope eps*mean(Y)+denorm_min. Rate uncertainty divides by this dt.
  void add_storage(std::size_t independent,double volume,double before_rho,double after_rho,
      long double before_y,long double after_y) noexcept {
    storage(mapping_[independent],volume,before_rho,after_rho,before_y,after_y,true);
  }
  void add_total_storage(double volume,double before_rho,double after_rho) noexcept {
    storage(rows_.size()-1,volume,before_rho,after_rho,1,1,false);
  }
  static bool admissible(const DriverCompositionBalance& row) noexcept {
    return std::isfinite(row.defect) && std::isfinite(row.relative_defect) &&
        std::isfinite(row.storage_roundoff_bound) && row.storage_roundoff_bound>=0 &&
        (row.relative_defect<1e-6 || std::abs(row.defect)<=row.storage_roundoff_bound);
  }
  void freeze_transport(std::size_t independent, long double value) noexcept {
    add(independent, transport, value / fields_);
    ++transport_counts_[mapping_[independent]];
  }
  void start_pressure() noexcept {
    for (auto& row : rows_) row[pressure] = 0;
    std::fill(pressure_counts_.begin(), pressure_counts_.end(), 0);
  }
  void add_pressure(std::size_t independent, long double value) noexcept {
    add(independent, pressure, value / fields_);
    ++pressure_counts_[mapping_[independent]];
  }
  static long double flux_sum(const CartesianKernelPlan& kernels,
      ConstFaceFluxView flux, Span<const std::uint8_t> activity) noexcept {
    const auto cells = kernels.cells();
    long double sum{}; std::size_t i{};
    for (int z=0; z<cells.z; ++z) for (int y=0; y<cells.y; ++y)
      for (int x=0; x<cells.x; ++x, ++i) {
        if (activity.size && !activity.data[i]) continue;
        const Int3 c{x,y,z};
        sum += static_cast<long double>(flux.x.unchecked({x+1,y,z})) - flux.x.unchecked(c);
        sum += static_cast<long double>(flux.y.unchecked({x,y+1,z})) - flux.y.unchecked(c);
        sum += static_cast<long double>(flux.z.unchecked({x,y,z+1})) - flux.z.unchecked(c);
      }
    return sum;
  }
  static long double pressure_sum(const CartesianKernelPlan& kernels,
      const BoundaryPlan& boundary, ConstFieldView q, ConstFieldView velocity,
      ConstFaceFluxView old_flux, ConstFaceFluxView flux,
      Span<const std::uint8_t> activity) noexcept {
    const auto cells = kernels.cells();
    const std::array<ConstFaceFieldView,3> a{old_flux.x,old_flux.y,old_flux.z};
    const std::array<ConstFaceFieldView,3> b{flux.x,flux.y,flux.z};
    long double sum{}; std::size_t i{};
    for (int z=0; z<cells.z; ++z) for (int y=0; y<cells.y; ++y)
      for (int x=0; x<cells.x; ++x, ++i) {
        if (activity.size && !activity.data[i]) continue;
        for (unsigned f=0; f<6; ++f) {
          const unsigned axis=f/2; Int3 face{x,y,z};
          if (f%2) ++(axis==0 ? face.x : axis==1 ? face.y : face.z);
          const double dm=b[axis].unchecked(face)-a[axis].unchecked(face);
          if (dm==0) continue;
          Int3 donor=face;
          if (dm>0) --(axis==0 ? donor.x : axis==1 ? donor.y : donor.z);
          sum += (f%2 ? 1.L : -1.L)*dm*scalar_upwind_donor(boundary,q,donor,velocity);
        }
      }
    return sum;
  }
  Status finish(ReductionEngine& reductions, DriverConservationReport& out) {
    Status local;
    if (!active()) local=invalid();
    for (const auto s : mapping_)
      if (transport_counts_[s]!=fields_ || pressure_counts_[s]!=fields_) local=invalid();
    for (const auto& row : rows_) for (const auto value : row)
      if (!std::isfinite(value)) local={StatusCode::numerical_failure,17861};
    for(const auto value:storage_bounds_)
      if(!std::isfinite(value) || value<0)local={StatusCode::numerical_failure,17861};
    auto status=reductions.consensus(local); if (!status) return status;
    auto& dependent=rows_[dependent_]; dependent=rows_.back();
    storage_bounds_[dependent_]=storage_bounds_.back();
    for (const auto s : mapping_) for (unsigned t=0;t<term_count;++t) dependent[t]-=rows_[s][t];
    for(const auto s:mapping_)storage_bounds_[dependent_]+=storage_bounds_[s];
    const auto ns=identity_->species_names.size(), ne=identity_->element_names.size();
    std::vector<Row> global;
    std::vector<DriverCompositionBalance> species, elements;
    std::vector<double> bounds;
    try {
      global.resize(ns);species.resize(ns);elements.resize(ne);
      bounds.resize(ns);
      for(std::size_t s=0;s<ns;++s)species[s].name=identity_->species_names[s];
      for(std::size_t e=0;e<ne;++e)elements[e].name=identity_->element_names[e];
    }
    catch (...) {local={StatusCode::allocation_failure,17861};}
    status=reductions.consensus(local); if (!status) return status;
    // Fill each reduction packet to the negotiated capacity. Species and
    // storage-roundoff bounds remain separate coordinates in the same SUM.
    const std::size_t width=term_count+1;
    std::vector<double> packed, reduced;
    try {
      packed.resize(ns*width);reduced.resize(ns*width);
      for(std::size_t s=0;s<ns;++s) {
        for(unsigned t=0;t<term_count;++t)packed[s*width+t]=static_cast<double>(rows_[s][t]);
        packed[s*width+term_count]=static_cast<double>(storage_bounds_[s]);
      }
    } catch(...) {local={StatusCode::allocation_failure,17861};}
    status=reductions.consensus(local);if(!status)return status;
    if(reductions.capacity()==0)return invalid();
    for(std::size_t first=0;first<packed.size();first+=reductions.capacity()) {
      const auto count=std::min(reductions.capacity(),packed.size()-first);
      status=reductions.checked_sum({packed.data()+first,count},{reduced.data()+first,count},{});
      if(!status)return status;
    }
    for(std::size_t s=0;s<ns;++s) {
      for(unsigned t=0;t<term_count;++t)global[s][t]=reduced[s*width+t];
      bounds[s]=reduced[s*width+term_count];
      report(species[s],global[s],bounds[s]);
    }
    for (std::size_t e=0;e<ne;++e) {
      Row element{};
      long double bound{};
      for (std::size_t s=0;s<ns;++s) {
        const long double weight=static_cast<long double>(identity_->element_counts[s*ne+e])/
            identity_->molecular_weights_kg_per_kmol[s];
        for (unsigned t=0;t<term_count;++t) element[t]+=weight*global[s][t];
        bound+=std::abs(weight)*bounds[s];
      }
      report(elements[e],element,bound);
    }
    out.composition_valid=true;out.composition_revision=generation_;
    out.composition_duration=dt_;out.composition_after_parcel_exchange=after_parcels_;
    out.species_balance=std::move(species);out.element_balance=std::move(elements);
    return {};
  }
 private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan,17861}; }
  static void report(DriverCompositionBalance& out,const Row& row,long double bound) noexcept {
    const long double defect=row[temporal]+row[transport]+row[pressure]-row[noise]-row[mixing]-row[chemistry];
    long double scale=1e-30L;
    for (unsigned t=temporal;t<term_count;++t) scale=std::max(scale,std::abs(row[t]));
    out.accepted_inventory=static_cast<double>(row[accepted]);
    out.current_inventory=static_cast<double>(row[current]);
    out.temporal_rate=static_cast<double>(row[temporal]);
    out.transport_outflow=static_cast<double>(row[transport]);
    out.pressure_outflow=static_cast<double>(row[pressure]);
    out.noise_source=static_cast<double>(row[noise]);
    out.mixing_source=static_cast<double>(row[mixing]);
    out.chemistry_source=static_cast<double>(row[chemistry]);
    out.defect=static_cast<double>(defect);
    out.relative_defect=static_cast<double>(std::abs(defect)/scale);
    out.storage_roundoff_bound=static_cast<double>(bound);
    out.roundoff_applied=out.relative_defect>=1e-6 && admissible(out);
  }
  void storage(std::size_t index,double volume,double before_rho,double after_rho,
      long double before_y,long double after_y,bool uncertain_y) noexcept {
    if(!std::isfinite(volume) || volume<=0 || !std::isfinite(before_rho) || before_rho<=0 ||
        !std::isfinite(after_rho) || after_rho<=0 || !std::isfinite(before_y) || before_y<0 ||
        !std::isfinite(after_y) || after_y<0) {
      storage_bounds_[index]=std::numeric_limits<long double>::quiet_NaN();return;
    }
    const auto spacing=[](long double value) {
      return std::numeric_limits<double>::epsilon()*std::abs(value)+
          std::numeric_limits<double>::denorm_min();
    };
    const auto bound=[&](long double rho,long double y) {
      const long double dv=spacing(volume),dr=spacing(rho),dy=uncertain_y ? spacing(y) : 0;
      return dv*(rho*y+rho*dy+dr*y+dr*dy)+
          static_cast<long double>(volume)*(rho*dy+dr*y+dr*dy);
    };
    const long double before=static_cast<long double>(volume)*before_rho*before_y;
    const long double after=static_cast<long double>(volume)*after_rho*after_y;
    auto& row=rows_[index];
    row[accepted]+=before;row[current]+=after;row[temporal]+=(after-before)/dt_;
    storage_bounds_[index]+=(bound(before_rho,before_y)+bound(after_rho,after_y))/dt_;
  }
  const portable::GasIdentity* identity_{};
  std::vector<Row> rows_;
  std::vector<long double> storage_bounds_;
  std::vector<std::size_t> mapping_,transport_counts_,pressure_counts_;
  std::size_t dependent_{},fields_{};
  RevisionToken generation_{}; double dt_{}; bool after_parcels_{};
};
} // namespace hundun::v04::detail

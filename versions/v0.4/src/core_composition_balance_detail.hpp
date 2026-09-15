// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "hundun/v04_app.hpp"
#include "hundun/v04_portable.hpp"
#include "solver_equation_detail.hpp"
#include "solver_scalar_boundary_detail.hpp"

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
    auto status=reductions.consensus(local); if (!status) return status;
    auto& dependent=rows_[dependent_]; dependent=rows_.back();
    for (const auto s : mapping_) for (unsigned t=0;t<term_count;++t) dependent[t]-=rows_[s][t];
    const auto ns=identity_->species_names.size(), ne=identity_->element_names.size();
    std::vector<Row> global;
    std::vector<DriverCompositionBalance> species, elements;
    try {
      global.resize(ns);species.resize(ns);elements.resize(ne);
      for(std::size_t s=0;s<ns;++s)species[s].name=identity_->species_names[s];
      for(std::size_t e=0;e<ne;++e)elements[e].name=identity_->element_names[e];
    }
    catch (...) {local={StatusCode::allocation_failure,17861};}
    status=reductions.consensus(local); if (!status) return status;
    for (std::size_t s=0;s<ns;++s) {
      std::array<double,term_count> input{}, result{};
      for (unsigned t=0;t<term_count;++t) input[t]=static_cast<double>(rows_[s][t]);
      status=reductions.checked_sum({input.data(),term_count},{result.data(),term_count},{});
      if (!status) return status;
      for (unsigned t=0;t<term_count;++t) global[s][t]=result[t];
      report(species[s],global[s]);
    }
    for (std::size_t e=0;e<ne;++e) {
      Row element{};
      for (std::size_t s=0;s<ns;++s) {
        const long double weight=static_cast<long double>(identity_->element_counts[s*ne+e])/
            identity_->molecular_weights_kg_per_kmol[s];
        for (unsigned t=0;t<term_count;++t) element[t]+=weight*global[s][t];
      }
      report(elements[e],element);
    }
    out.composition_valid=true;out.composition_revision=generation_;
    out.composition_duration=dt_;out.composition_after_parcel_exchange=after_parcels_;
    out.species_balance=std::move(species);out.element_balance=std::move(elements);
    return {};
  }
 private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan,17861}; }
  static void report(DriverCompositionBalance& out,const Row& row) noexcept {
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
  }
  const portable::GasIdentity* identity_{};
  std::vector<Row> rows_;
  std::vector<std::size_t> mapping_,transport_counts_,pressure_counts_;
  std::size_t dependent_{},fields_{};
  RevisionToken generation_{}; double dt_{}; bool after_parcels_{};
};
} // namespace hundun::v04::detail

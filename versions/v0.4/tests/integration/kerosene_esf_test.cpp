// SPDX-License-Identifier: Apache-2.0
// Small native 624CF gas fixture: full histories and absolute face-flux audit.
#include "hundun/v04_app.hpp"
#include "../../src/models_tcr_dynamic_detail.hpp"
#include "../../src/core_tcr_dyn711_history_detail.hpp"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string_view>
#include <iostream>
using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);struct End{~End(){MPI_Finalize();}} end;int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const auto all=[&](bool value){int ok=value;MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);return bool(ok);};
  if(argc!=4 && (argc!=5 || (std::string_view(argv[4])!="--wall" && std::string_view(argv[4])!="--evaporated" && std::string_view(argv[4])!="--collapsed" && std::string_view(argv[4])!="--reference-linear")))return 2;
  const bool collapsed=argc==5 && std::string_view(argv[4])=="--collapsed";
  const bool reference=argc==5 && std::string_view(argv[4])=="--reference-linear";
  const double field_limit=reference ? 1e-4 : 1e-11;
  const bool wall=argc==5 && !collapsed && !reference;
  const bool evaporated=wall && std::string_view(argv[4])=="--evaporated";
  ValidatedModel model;auto s=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],model);
  CompiledCasePlan plan;if(s)s=ProductCompiler::compile(MPI_COMM_WORLD,model,argv[1],plan);
  ProductDriver driver;if(s)s=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;if(s)s=driver.restart_expected(expected);
  RestartImage a,b;if(s)s=RestartReader::load(MPI_COMM_WORLD,argv[2],expected,a);
  if(s)s=RestartReader::load(MPI_COMM_WORLD,argv[3],expected,b);
  if(!all(bool(s))) {if(!rank)std::cout<<"load "<<unsigned(s.code)<<'/'<<s.detail<<std::endl;return 3;}
  // Parcel ownership can leave a valid reader rank with an empty payload.
  int global_records = !a.cell_records.empty();
  MPI_Allreduce(MPI_IN_PLACE,&global_records,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  bool records_equal=a.cell_records==b.cell_records;
  double history_error{}, kappa_roundoff_bound{};
  if(model.reaction.esf && model.reaction.esf->tcr.model==TcrModel::cdphyso_dynamic_v1 &&
      a.cell_record_bytes==b.cell_record_bytes && a.cell_records.size()==b.cell_records.size() &&
      a.cell_record_lengths==b.cell_record_lengths) {
    const auto real=[](const std::uint8_t* p) {
      std::uint64_t bits{};for(unsigned i=0;i<8;++i)bits|=std::uint64_t(p[i])<<(8*i);
      double v;std::memcpy(&v,&bits,8);return v;
    };
    const auto ns=model.thermophysics.species.size();
    const auto prefix=model.spray ? 24U : 0U;
    const auto tcr_width=24+8*(5*ns+3);
    const auto count=model.spray ? a.cell_record_lengths.size() :
        (a.cell_record_bytes ? a.cell_records.size()/a.cell_record_bytes : 0);
    records_equal=count>0;std::size_t offset{};
    for(std::size_t cell=0;cell<count;++cell) {
      const auto width=model.spray ? a.cell_record_lengths[cell] : a.cell_record_bytes;
      if(width<prefix+tcr_width || offset>a.cell_records.size() || width>a.cell_records.size()-offset) {records_equal=false;break;}
      const auto *ap=a.cell_records.data()+offset,*bp=b.cell_records.data()+offset;
      records_equal &= std::equal(ap,ap+prefix+24,bp) &&
          std::equal(ap+prefix+tcr_width,ap+width,bp+prefix+tcr_width);
      for(std::size_t j=0;j<5*ns+3;++j) {
        const auto *av=ap+prefix+24+8*j,*bv=bp+prefix+24+8*j;
        if(j<5*ns && j%5==4)records_equal &= std::equal(av,av+8,bv);
        else {
          const double x=real(av),y=real(bv);
          const double gap=std::abs(x-y)/std::max({1.,std::abs(x),std::abs(y)});
          double tolerance=1e-11;
          if(j<5*ns && (j%5==2 || j%5==3)) {
            const auto species=j/5;
            const auto base=prefix+24+8*(5*species);
            const double pdf=real(ap+base), psr=real(ap+base+8);
            const auto select=[&](double f,double s) {
              const auto c=tcr::detail::cdphyso_species_control(.3,f,s,
                  model.reaction.esf->tcr.weak_rate_threshold);
              return j%5==2 ? c.selected : c.effective;
            };
            // Rates subtract two FP64 mass fractions before division by dt*W.
            // Bound their cancellation uncertainty before the nonlinear root;
            // raw rates, Cd, fields and every discrete state retain 1e-11/exact checks.
            const double ulp_rate=4*std::numeric_limits<double>::epsilon()/
                (a.dt*model.thermophysics.species[species].molecular_weight);
            const auto own=tcr::detail::cdphyso_species_control(.3,pdf,psr,
                model.reaction.esf->tcr.weak_rate_threshold);
            if(own.state==tcr::detail::SpeciesControlState::direct ||
               own.state==tcr::detail::SpeciesControlState::projected) {
              const double center=select(pdf,psr);
              records_equal &= std::abs(x-center)<=1e-11;
              for(int df:{-1,0,1})for(int ds:{-1,0,1})
                tolerance=std::max(tolerance,std::abs(select(pdf+df*ulp_rate,psr+ds*ulp_rate)-center)+1e-11);
              const double other=select(real(bp+base),real(bp+base+8));
              records_equal &= std::abs(y-other)<=1e-11;
              kappa_roundoff_bound=std::max(kappa_roundoff_bound,tolerance);
            }
          }
          if(gap>=tolerance && gap>history_error)
            std::printf("tcr_history_difference rank=%d cell=%zu coordinate=%zu first=%.17g second=%.17g tolerance=%.17g\n",rank,cell,j,x,y,tolerance);
          records_equal &= std::isfinite(x) && std::isfinite(y) && gap<tolerance;
          history_error=std::max(history_error,gap);
        }
      }
      offset+=width;
    }
    records_equal &= offset==a.cell_records.size();
  }
  if(model.reaction.esf && model.reaction.esf->tcr.model==TcrModel::dyn711_v1 &&
      a.cell_record_bytes==b.cell_record_bytes && a.cell_records.size()==b.cell_records.size()) {
    const auto ns=model.thermophysics.species.size(),width=40+40*ns;
    const auto real=[](const std::uint8_t *p) {
      std::uint64_t bits{};for(unsigned j=0;j<8;++j)bits|=std::uint64_t(p[j])<<(8*j);
      double value;std::memcpy(&value,&bits,8);return value;
    };
    records_equal=!model.spray && a.cell_record_bytes==width && a.cell_records.size()%width==0;
    if(records_equal) {
      detail::Dyn711History left,right;
      const auto count=a.cell_records.size()/width;
      left.configure(1,count,ns,2);right.configure(1,count,ns,2);
      records_equal=bool(left.restore({a.cell_records.data(),a.cell_records.size()},a.step)) &&
          bool(right.restore({b.cell_records.data(),b.cell_records.size()},b.step));
      left.commit();right.commit();
      records_equal &= left.statistics_calls()==a.step && right.statistics_calls()==b.step;
      for(std::size_t i=0;i<count;++i) {
        const auto *x=a.cell_records.data()+i*width,*y=b.cell_records.data()+i*width;
        records_equal &= std::equal(x,x+32,y);
        for(std::size_t j=32;j<width;j+=8) {
          if(j>=40 && (j-40)%40==32)records_equal &= std::equal(x+j,x+j+8,y+j);
          else {
            const double av=real(x+j),bv=real(y+j);
            const double error=std::abs(av-bv)/std::max({1.,std::abs(av),std::abs(bv)});
            records_equal &= std::isfinite(av) && std::isfinite(bv) && error<1e-11;
            if(error>1e-11 && error>history_error)
              std::printf("dyn711_history_difference rank=%d cell=%zu byte=%zu left=%.17g right=%.17g relative=%.17g global=%zu,%zu,%zu\n",
                  rank,i,j,av,bv,error,i%a.patch.cells.x+a.patch.begin.x,
                  (i/a.patch.cells.x)%a.patch.cells.y+a.patch.begin.y,
                  i/(a.patch.cells.x*a.patch.cells.y)+a.patch.begin.z);
            history_error=std::max(history_error,error);
          }
        }
      }
    }
    if(!rank)std::printf("dyn711_history clocks=restored rates_roots_cphi_limit=1e-11 valid=%d\n",int(records_equal));
  }
  MPI_Allreduce(MPI_IN_PLACE,&history_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&kappa_roundoff_bound,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(!rank && model.reaction.esf && dynamic_tcr_model(model.reaction.esf->tcr.model))
    std::printf("dynamic_tcr_compare discrete=exact parcels=exact history_relative=%.17g rate_cd_limit=1e-11 kappa_roundoff_bound=%.17g\n",history_error,kappa_roundoff_bound);
  const bool metadata=a.step==b.step && a.time==b.time && a.dt==b.dt && a.method_history_signature==b.method_history_signature && a.plan==b.plan && a.schema==b.schema && a.geometry==b.geometry && a.source_format_version==b.source_format_version && !a.backward_euler_recovery && !b.backward_euler_recovery && a.controller_state==b.controller_state && a.pressure_reference==b.pressure_reference && a.previous_pressure_reference==b.previous_pressure_reference && a.closed_mass_target==b.closed_mass_target && records_equal &&
      a.cell_record_lengths==b.cell_record_lengths &&
      a.cell_record_identity==b.cell_record_identity && a.cell_record_bytes==b.cell_record_bytes &&
      (!model.spray || (global_records && a.cell_record_identity!=0));
  if(!all(metadata)) {
    if(!metadata)std::cout << "metadata mismatch rank=" << rank << " step=" << (a.step==b.step)
      << " time=" << (a.time==b.time) << " dt=" << (a.dt==b.dt)
      << " method=" << (a.method_history_signature==b.method_history_signature)
      << " plan=" << (a.plan==b.plan) << " schema=" << (a.schema==b.schema)
      << " geometry=" << (a.geometry==b.geometry)
      << " source=" << (a.source_format_version==b.source_format_version)
      << " recovery=" << (!a.backward_euler_recovery && !b.backward_euler_recovery)
      << " controller=" << (a.controller_state==b.controller_state)
      << " pressure=" << (a.pressure_reference==b.pressure_reference && a.previous_pressure_reference==b.previous_pressure_reference)
      << " mass=" << (a.closed_mass_target==b.closed_mass_target)
      << " records=" << (a.cell_records==b.cell_records)
      << " lengths=" << (a.cell_record_lengths==b.cell_record_lengths)
      << " record_identity=" << (a.cell_record_identity==b.cell_record_identity)
      << " record_bytes=" << (a.cell_record_bytes==b.cell_record_bytes)
      << " spray=" << bool(model.spray) << " stored_bytes=" << a.cell_records.size()
      << " identity=" << a.cell_record_identity << std::endl;
    // Continue with field diagnostics; metadata still participates in the final result.
  }
  if(model.reaction.esf->tcr.model==TcrModel::cdphyso_dynamic_v1) {
    const auto integer=[](const std::uint8_t* p,unsigned width) {
      std::uint64_t v{};for(unsigned i=0;i<width;++i)v|=std::uint64_t(p[i])<<(8*i);return v;
    };
    const auto real=[&](const std::uint8_t* p) {
      const auto bits=integer(p,8);double v;std::memcpy(&v,&bits,8);return v;
    };
    const auto ns=model.thermophysics.species.size();
    std::size_t offset{};bool valid=true;double live[3]{};
    const auto prefix=model.spray ? 24U : 0U;
    const auto count=model.spray ? a.cell_record_lengths.size() :
        (a.cell_record_bytes ? a.cell_records.size()/a.cell_record_bytes : 0);
    valid &= count>0;
    for(std::size_t cell=0;cell<count;++cell) {
      const auto width=model.spray ? a.cell_record_lengths[cell] : a.cell_record_bytes;
      if(width<prefix+24+8*(5*ns+3) || offset>a.cell_records.size() || width>a.cell_records.size()-offset) {valid=false;break;}
      const auto *p=a.cell_records.data()+offset+prefix;
      valid &= integer(p,8)==a.step && integer(p+8,8)==a.step%4 && integer(p+16,8)==ns;
      for(std::size_t q=0;q<ns;++q) {
        live[0]=std::max(live[0],std::abs(real(p+24+8*(5*q))));
        live[1]=std::max(live[1],1-real(p+24+8*(5*q+3)));
      }
      for(unsigned g=0;g<3;++g)live[2]=std::max(live[2],
          std::abs(real(p+24+8*(5*ns+g))-2/model.reaction.mixing_c_z));
      offset+=width;
    }
    MPI_Allreduce(MPI_IN_PLACE,live,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(!rank)std::printf("dynamic_tcr history=audited cadence_phase=%llu rate_max=%.17g kappa_change=%.17g cd_change=%.17g\n",
        static_cast<unsigned long long>(a.step%4),live[0],live[1],live[2]);
    if(!all(valid && live[0]>0 && live[1]>0 && live[2]>0))return 10;
  }
  if(wall) {
    // Native V5 variable cell records: header, TCR, parcels, injectors.
    // Check the actual trajectory after the specified cube-wall encounter.
    const auto integer=[](const std::uint8_t* p,unsigned width) {
      std::uint64_t value{};for(unsigned i=0;i<width;++i)value|=std::uint64_t(p[i])<<(8*i);
      return value;
    };
    const auto real=[&](const std::uint8_t* p) {
      const auto bits=integer(p,8);double value;std::memcpy(&value,&bits,8);return value;
    };
    bool valid=model.spray.has_value() && model.immersed_boundary.has_value() && a.source_format_version==5;
    std::size_t offset{};unsigned long long parcels{};
    for(const auto bytes:a.cell_record_lengths) {
      if(!bytes)continue;
      if(bytes<24 || offset>a.cell_records.size() || bytes>a.cell_records.size()-offset) {valid=false;break;}
      const auto* record=a.cell_records.data()+offset;
      const auto count=integer(record+8,4),injectors=integer(record+12,4),tcr=integer(record+16,4);
      const auto version=integer(record+20,4);
      const unsigned width=version==2 ? 192 : 144;
      if((version!=1 && version!=2) || 24+tcr+width*count+24*injectors!=bytes) {valid=false;break;}
      const auto* p=record+24+tcr;
      for(std::uint64_t i=0;i<count;++i,p+=width) {
        const double x=real(p+16),u=real(p+40);
        const double length=model.mesh.upper.x-model.mesh.lower.x;
        const double wall_x=model.mesh.lower.x+.375*length;
        const double travel=200.*a.dt*a.step+1e-6*length;
        if(!(x<wall_x && x>wall_x-travel && std::isfinite(u) && u<0))
          std::printf("wall_probe rank=%d x=%.17g u=%.17g wall=%.17g travel=%.17g\n",rank,x,u,wall_x,travel);
        valid &= x<wall_x && x>wall_x-travel && std::isfinite(u) && u<0;
        if(model.spray->breakup==SprayBreakupModel::stochastic_sgs) {
          // Exact cross-partition record comparison above includes the full
          // lineage; require live exposure history in this native run too.
          valid &= version==2 && integer(p+144,8)==1 &&
                   real(p+160)>0 && std::isfinite(real(p+152));
        }
        ++parcels;
      }
      offset+=bytes;
    }
    valid &= offset==a.cell_records.size();
    MPI_Allreduce(MPI_IN_PLACE,&parcels,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    if(!all(valid && parcels==(evaporated ? 0 : a.step)))return 9;
    if(!rank)std::printf("ibm_parcel_wall retained=%llu complete_evaporation=%d records=exact\n",parcels,int(evaporated));
  }
  bool passed=metadata;double worst{};
  const auto role_count=[&](RestartFieldRole role) {
    return std::count_if(a.fields.begin(),a.fields.end(),[&](const auto& f){return f.role==role;});
  };
  if(!all(role_count(RestartFieldRole::stochastic_field)==model.reaction.esf->fields &&
          role_count(RestartFieldRole::stochastic_auxiliary)==1))return 8;
  if(collapsed) {
    const RestartImageField *first=nullptr;
    bool same=true;double spread{};
    for(const auto &field:a.fields)if(field.role==RestartFieldRole::stochastic_field) {
      if(!first) {first=&field;continue;}
      if(field.values.size()!=first->values.size()) {same=false;continue;}
      for(std::size_t i=0;i<field.values.size();++i) {
        const double x=field.values[i],y=first->values[i];
        const double error=std::abs(x-y)/std::max({1.,std::abs(x),std::abs(y)});
        same &= std::isfinite(x) && std::isfinite(y) && error<1e-13;
        spread=std::max(spread,error);
      }
    }
    const bool valid=all(same && first && model.reaction.esf->tcr.model==TcrModel::dyn711_v1);
    MPI_Allreduce(MPI_IN_PLACE,&spread,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(!rank)std::printf("dyn711_remix physical_fields_collapsed=%d relative_spread=%.17g\n",int(valid),spread);
    if(!valid)return 10;
  }
  // Mechanical pressure is a perturbation of the stored absolute reference.
  // Its spatial work U.grad(p) propagates primary pressure differences into
  // the persisted enthalpy rate. Bound that propagation from the actual
  // compared states and mesh, independently of the measured rate difference.
  const auto pressure_work_bound=[&](const std::vector<RestartImageField>& left,
                                      const std::vector<RestartImageField>& right) {
    double gap{}, speed{};
    for(std::size_t f=0;f<left.size();++f) {
      if(f>=right.size() || left[f].values.size()!=right[f].values.size())continue;
      if(left[f].role==RestartFieldRole::pressure_perturbation)
        for(std::size_t i=0;i<left[f].values.size();++i)
          gap=std::max(gap,std::abs(left[f].values[i]-right[f].values[i]));
      if(left[f].role==RestartFieldRole::velocity)
        for(std::size_t i=0;i+2<left[f].values.size();i+=3)
          speed=std::max({speed,std::abs(left[f].values[i])+std::abs(left[f].values[i+1])+std::abs(left[f].values[i+2]),
              std::abs(right[f].values[i])+std::abs(right[f].values[i+1])+std::abs(right[f].values[i+2])});
    }
    double values[]{gap,speed};MPI_Allreduce(MPI_IN_PLACE,values,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    double dx=INFINITY;
    const double lo[]{model.mesh.lower.x,model.mesh.lower.y,model.mesh.lower.z};
    const double hi[]{model.mesh.upper.x,model.mesh.upper.y,model.mesh.upper.z};
    const int n[]{a.global_cells.x,a.global_cells.y,a.global_cells.z};
    for(unsigned d=0;d<3;++d) {
      const auto& faces=model.mesh.runtime_faces[d];
      if(faces.empty())dx=std::min(dx,(hi[d]-lo[d])/n[d]);
      else for(std::size_t i=1;i<faces.size();++i)dx=std::min(dx,faces[i]-faces[i-1]);
    }
    return 2*values[0]*values[1]/dx;
  };
  const double current_work_bound=pressure_work_bound(a.fields,b.fields);
  const double previous_work_bound=pressure_work_bound(a.previous_fields,b.previous_fields);
  double maximum_flux_difference{};unsigned long long total{};
  const auto fields=[&](const char* level,const std::vector<RestartImageField>& x,const std::vector<RestartImageField>& y) {
    if(!all(x.size()==y.size()))return false;
    for(std::size_t f=0;f<x.size();++f) {
      if(!all(x[f].field==y[f].field && x[f].role==y[f].role && x[f].components==y[f].components && x[f].values.size()==y[f].values.size()))return false;
      for(unsigned c=0;c<x[f].components;++c) {
        double v[2]{};unsigned long long differences{};
        for(std::size_t i=c;i<x[f].values.size();i+=x[f].components) {
          const double aa=x[f].values[i],bb=y[f].values[i];v[0]=std::max(v[0],std::abs(aa-bb));v[1]=std::max({v[1],std::abs(aa),std::abs(bb)});differences+=aa!=bb;
          if(!std::isfinite(aa)||!std::isfinite(bb))v[0]=INFINITY;
          ++total;
        }
        MPI_Allreduce(MPI_IN_PLACE,v,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&differences,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
        double scale=std::max(1.,v[1]);
        if(x[f].role==RestartFieldRole::pressure_perturbation)
          scale=std::max({scale,std::abs(a.pressure_reference),std::abs(a.previous_pressure_reference)});
        const double propagated=x[f].role==RestartFieldRole::enthalpy_nonadvective_rate ?
            (std::string_view(level)=="rate" ? current_work_bound : previous_work_bound) : 0.;
        const double relative=v[0]/scale;worst=std::max(worst,relative);
        passed &= v[0]<field_limit*scale+propagated;
        if(!rank && propagated>0)std::printf("pressure_work_propagation level=%s absolute_bound=%.17g\n",level,propagated);
        if(!rank)std::printf("compare level=%s field=%u role=%u component=%u abs=%.17g scale=%.17g relative=%.17g unequal=%llu\n",level,x[f].field,unsigned(x[f].role),c,v[0],scale,relative,differences);
      }
    }
    return true;
  };
  if(!fields("current",a.fields,b.fields) || !fields("previous",a.previous_fields,b.previous_fields) || !fields("rate",a.accepted_rate_fields,b.accepted_rate_fields) || !fields("previous_rate",a.previous_rate_fields,b.previous_rate_fields))return 5;
  double flux_scale=0.;
  for(unsigned axis=0;axis<3;++axis)for(double value:a.final_mass_flux[axis])flux_scale=std::max(flux_scale,std::abs(value));
  MPI_Allreduce(MPI_IN_PLACE,&flux_scale,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  const double flux_limit=reference ? 1e-4*std::max(1e-30,flux_scale) : 1e-11;
  for(unsigned level=0;level<2;++level)for(unsigned axis=0;axis<3;++axis) {
    const auto& x=level ? a.previous_mass_flux[axis] : a.final_mass_flux[axis];const auto& y=level ? b.previous_mass_flux[axis] : b.final_mass_flux[axis];
    if(!all(x.size()==y.size()))return 6;
    double v[2]{};
    for(std::size_t i=0;i<x.size();++i){v[0]=std::max(v[0],std::abs(x[i]-y[i]));v[1]=std::max({v[1],std::abs(x[i]),std::abs(y[i])});if(!std::isfinite(x[i]) || !std::isfinite(y[i]))v[0]=INFINITY;++total;}
    MPI_Allreduce(MPI_IN_PLACE,v,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    const double relative=v[0]/std::max(1e-30,v[1]);maximum_flux_difference=std::max(maximum_flux_difference,v[0]);passed &= v[0]<flux_limit;
    if(!rank)std::printf("compare_flux level=%u axis=%u abs=%.17g scale=%.17g relative=%.17g\n",level,axis,v[0],v[1],relative);
  }
  MPI_Allreduce(MPI_IN_PLACE,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if(!rank)std::printf("kerosene_esf_compare step=%llu dt=%.17g metadata=exact field_worst=%.17g field_limit=%.17g flux_abs_max_kg_s=%.17g flux_abs_limit_kg_s=%.17g values=%llu passed=%d\n",static_cast<unsigned long long>(a.step),a.dt,worst,field_limit,maximum_flux_difference,flux_limit,total,int(passed));
  return passed ? 0 : 7;
}

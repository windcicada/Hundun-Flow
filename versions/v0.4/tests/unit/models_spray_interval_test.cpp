// SPDX-License-Identifier: Apache-2.0
// THICK_EX interval, actual inventory exchange and adaptive event checks.
#include "models_spray_events_detail.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace hundun::v04;
using namespace spray;
using namespace spray::detail;
namespace {
struct Geometry final : ParcelEventGeometryProvider {
  ParcelEventQueryReport query(const SprayParcelState &, const SprayParcelState &,
      double, double, ParcelPass, portable::Revision revision) const noexcept override {
    ParcelEventQueryReport out; out.available=true; out.revision=revision; return out;
  }
};
struct IncrementError final : ParcelIntervalProvider {
  const ParcelIntervalProvider &interval;
  explicit IncrementError(const ParcelIntervalProvider &value) : interval(value) {}
  ParcelIntervalReport advance(const SprayParcelState &p, double t, double dt,
      ParcelPass pass, portable::Revision revision) const noexcept override {
    auto out=interval.advance(p,t,dt,pass,revision);
    out.exchange_from_endpoints=false;
    return out;
  }
};
struct Environment final : ParcelTransferEnvironmentProvider {
  LiquidPropertyPack pack{};
  LiquidPropertyService liquid{&pack, 1};
  double far_temperature{350};
  bool stale{};
  Vector3 gas_velocity{};
  Environment() {
    pack.material_fingerprint = 91;
    pack.minimum_temperature_k = 200;
    pack.maximum_temperature_k = 650;
    pack.density_kg_per_m3.c[0] = 750;
    pack.density_kg_per_m3.kind = TemperatureCorrelationKind::polynomial_cubic;
    pack.density_kg_per_m3.reference_temperature_k = 350;
    pack.density_kg_per_m3.c[1] = -.2;
    pack.cp_j_per_kg_k.c[0] = 2200;
    pack.cp_j_per_kg_k.kind = TemperatureCorrelationKind::polynomial_cubic;
    pack.cp_j_per_kg_k.reference_temperature_k = 350;
    pack.cp_j_per_kg_k.c[1] = 1.5;
    pack.latent_heat_j_per_kg.c[0] = 250000;
    pack.surface_tension_n_per_m.c[0] = .025;
    pack.viscosity_pa_s.c[0] = .0008;
    pack.saturation_pressure.antoine_a = 4;
    pack.saturation_pressure.pressure_scale_pa = 1;
  }
  ParcelTransferEnvironment sample(const SprayParcelState &s, double,
      ParcelPass, portable::Revision revision) const noexcept override {
    ParcelTransferEnvironment out;
    out.available = true;
    out.revision = revision;
    if (stale) ++out.revision.input_revision;
    out.far_gas_dynamic_viscosity_pa_s = 2e-5;
    out.boiling_temperature_k = 500;
    out.vapor_prandtl_number = .7;
    out.liquid = &liquid;
    const double change = s.temperature_k - 350;
    out.liquid_absolute_enthalpy_j_per_kg = 770000 + 2200*change + .75*change*change;
    out.far_gas_density_kg_per_m3 = 1;
    auto &g = out.gas;
    g.gas_velocity_m_per_s = gas_velocity;
    g.gas_temperature_k = far_temperature;
    g.gas_vapor_mass_fraction = .001;
    g.thermodynamic_pressure_pa = 1e5;
    g.film_density_kg_per_m3 = 1;
    g.film_dynamic_viscosity_pa_s = 2e-5;
    g.film_thermal_conductivity_w_per_m_k = .05;
    g.film_vapor_diffusivity_m2_per_s = 2e-5;
    g.film_cp_j_per_kg_k = 1100;
    g.vapor_molecular_weight_kg_per_kmol = 28;
    g.carrier_molecular_weight_kg_per_kmol = 28;
    g.vapor_absolute_thermochemical_enthalpy_j_per_kg =
        1020000 + 2200*change;
    return out;
  }
};
}
int main() {
  Environment gas;
  FixedThickParcelIntervalProvider production_interval(gas);
  SprayParcelState state;
  state.id = {1,2};
  state.droplet_diameter_m = 1e-4;
  state.droplet_mass_kg = std::acos(-1.)/6*750*1e-12;
  state.temperature_k = 350;
  state.multiplicity = 3;
  state.liquid_material_fingerprint = 91;
  const portable::Revision revision{1,1,1};
  bool passed = true;
  for (double temperature : {350.,900.}) {
    gas.far_temperature = temperature;
    const double dt = temperature==350. ? 1e-6 : 1e-9;
    auto actual = production_interval.advance(state,0,dt,ParcelPass::corrector,revision);
    ThickExchangeInput input;
    input.parcel = state;
    input.liquid_properties = gas.liquid.evaluate({91,350}).properties;
    input.gas_temperature_k = temperature;
    input.gas_cp_j_per_kg_k = 1100;
    input.gas_dynamic_viscosity_pa_s = 2e-5;
    input.prandtl_number = .7;
    input.boiling_temperature_k = 500;
    input.vapor_absolute_thermochemical_enthalpy_j_per_kg = 1020000;
    input.duration_s = dt;
    auto expected = evaluate_thick_exchange(input);
    if (!actual.available || !expected.succeeded()) return 2;
    const double mass_error = std::abs(actual.parcel.droplet_mass_kg -
        expected.candidate_parcel.droplet_mass_kg)/state.droplet_mass_kg;
    const double temperature_error = std::abs(actual.parcel.temperature_k -
        expected.candidate_parcel.temperature_k);
    std::cout << std::setprecision(17) << "gas_T=" << temperature
        << " production_mass=" << actual.parcel.droplet_mass_kg
        << " THICK_EX_mass=" << expected.candidate_parcel.droplet_mass_kg
        << " relative_mass_difference=" << mass_error
        << " temperature_difference_K=" << temperature_error << '\n';
    passed &= mass_error <= 1e-11 && temperature_error <= 1e-8;
  }
  // At equilibrium the Stokes limit has an independent BE velocity solution.
  gas.far_temperature = state.temperature_k;
  gas.gas_velocity = {1e-8, -2e-8, 3e-8};
  const double dt = 1e-4;
  const double response = 18*2e-5*dt/(750*state.droplet_diameter_m*state.droplet_diameter_m);
  auto r = production_interval.advance(state,0,dt,ParcelPass::corrector,revision);
  if (!r.available) return 3;
  for (unsigned d=0;d<3;++d) {
    // Re ~ 1e-7; the retained SN correction is separately resolved below.
    const auto drag = evaluate_schiller_naumann_drag({gas.gas_velocity,{},1,2e-5,
        state.droplet_diameter_m,state.droplet_mass_kg,state.multiplicity,0});
    const double correction = drag.acceleration_m_per_s2[d]*dt/gas.gas_velocity[d];
    const double velocity = gas.gas_velocity[d]*correction/(1+correction);
    passed &= std::abs(r.parcel.velocity_m_per_s[d]-velocity)<1e-22;
    passed &= std::abs(correction/response-1)<1e-5;
    passed &= r.exchange.gas_momentum_delta_kg_m_per_s[d] ==
        -r.exchange.parcel_momentum_delta_kg_m_per_s[d];
  }
  // Absolute liquid h differs from the frozen gas h_v-L; exchange follows
  // actual inventory and remains invariant under consistent mass weighting.
  gas.gas_velocity = {2,-3,4}; gas.far_temperature = 900;
  r = production_interval.advance(state,0,1e-9,ParcelPass::corrector,revision);
  if (!r.available) return 4;
  const double change=r.parcel.temperature_k-350;
  const long double h0=770000, h1=770000+2200*change+.75*change*change;
  const double dh=static_cast<double>(state.multiplicity*
      ((static_cast<long double>(r.parcel.droplet_mass_kg)-state.droplet_mass_kg)*h0+
       r.parcel.droplet_mass_kg*(h1-h0)));
  passed &= std::abs(r.exchange.thermal_exchange_to_gas_j+dh)<1e-16;
  passed &= r.exchange.thermal_exchange_state_residual_j==0 &&
      r.exchange.mass_closure_residual_kg==0;
  const double density=750-.2*change;
  passed &= std::abs(density*std::acos(-1.)/6*std::pow(r.parcel.droplet_diameter_m,3)/
                    r.parcel.droplet_mass_kg-1)<1e-13;
  auto single=state;single.multiplicity=1;
  const auto one=production_interval.advance(single,0,1e-9,ParcelPass::corrector,revision);
  passed &= one.available && one.parcel.temperature_k==r.parcel.temperature_k &&
      std::abs(one.exchange.gas_mass_delta_kg*3-r.exchange.gas_mass_delta_kg)<1e-25;
  // Exact terminal event deposits the remaining inventory once.
  const auto terminal=production_interval.advance(state,0,1e-6,ParcelPass::corrector,revision);
  passed &= terminal.available && terminal.complete_evaporation &&
      terminal.elapsed_duration_s<1e-6 && terminal.parcel.droplet_mass_kg==0 &&
      terminal.parcel.droplet_diameter_m==0 &&
      terminal.exchange.gas_mass_delta_kg==state.droplet_mass_kg*state.multiplicity;
  ParcelEventsInput input;
  Geometry geometry;input.geometry=&geometry;
  input.accepted_parcel=state;input.revision=revision;input.interval=&production_interval;
  input.duration_s=1e-9;input.initial_substep_s=1e-9;input.minimum_substep_s=1e-14;
  const auto adaptive=integrate_parcel_events(input);
  passed &= adaptive.available && adaptive.parcel.droplet_mass_kg<state.droplet_mass_kg &&
      adaptive.exchange.thermal_exchange_state_residual_j==0;
  IncrementError increment(production_interval);input.interval=&increment;
  const auto redundant=integrate_parcel_events(input);
  passed &= !redundant.available && redundant.status==ParcelEventsStatus::capacity_exceeded;
  // Time refinement of the original BE phase relation. A dense trajectory
  // supplies the limit; its operator is separately compared with full Fortran.
  const auto trajectory = [&](unsigned n) {
    auto p=state;
    for (unsigned i=0;i<n;++i) {
      const auto next=production_interval.advance(p,i*1e-9/n,1e-9/n,
                                                  ParcelPass::corrector,revision);
      if (!next.available) {p.droplet_mass_kg=0;break;}
      p=next.parcel;
    }
    return p;
  };
  const auto dense=trajectory(8192), coarse=trajectory(32), fine=trajectory(64);
  const double error32=std::abs(coarse.droplet_mass_kg-dense.droplet_mass_kg)/state.droplet_mass_kg;
  const double error64=std::abs(fine.droplet_mass_kg-dense.droplet_mass_kg)/state.droplet_mass_kg;
  const double order=std::log(error32/error64)/std::log(2.);
  passed &= dense.droplet_mass_kg>0 && order>.9 && order<1.1;
  const double adaptive_error=std::abs(adaptive.parcel.droplet_mass_kg-dense.droplet_mass_kg)/state.droplet_mass_kg;
  passed &= adaptive_error<1e-4;
  std::cout<<"BE_mass_order="<<order<<" error32="<<error32<<" error64="<<error64
      <<" adaptive_mass_error="<<adaptive_error<<" segments="<<adaptive.segment_count
      <<" prior_failure_segment="<<redundant.failure_segment<<'\n';
  gas.stale=true;
  passed &= !production_interval.advance(state,0,1e-9,ParcelPass::corrector,revision).available;
  gas.stale=false;
  passed &= !production_interval.advance(state,0,-1,ParcelPass::corrector,revision).available;
  std::cout<<"thick_interval_inventory_drag_event passed="<<passed
      <<" adaptive="<<adaptive.available<<" status="<<int(adaptive.status)<<'\n';
  return passed ? 0 : 1;
}

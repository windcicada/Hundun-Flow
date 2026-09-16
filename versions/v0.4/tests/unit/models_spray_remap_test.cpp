// SPDX-License-Identifier: Apache-2.0
#include "models_spray_remap_detail.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
using namespace hundun::v04::spray::detail;
namespace {
void require(bool ok, const char* message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}
LiquidAsset material() {
  LiquidAsset a;
  a.pack.material_fingerprint = 17;
  a.pack.minimum_temperature_k = 200;
  a.pack.maximum_temperature_k = 680;
  a.reference_temperature_k = 300;
  a.reference_liquid_enthalpy_j_per_kg = -1e6;
  a.pack.density_kg_per_m3.c[0] = 800;
  a.pack.cp_j_per_kg_k.c[0] = 2000;
  a.pack.latent_heat_j_per_kg.c[0] = 3e5;
  a.pack.surface_tension_n_per_m.c[0] = .025;
  a.pack.viscosity_pa_s.c[0] = .001;
  a.pack.saturation_pressure.pressure_scale_pa = 1000;
  return a;
}
}
int main() {
  auto a = material();
  constexpr double m = 8e-12;
  for (double t : {200.,201.,300.,417.5,679.,680.}) {
    const double h = -1e6+2000*(t-300);
    const auto r = remap_liquid_mass_enthalpy(a,m,h);
    require(r.available && std::abs(r.temperature_k-t)<1e-11,
            "analytic constant-cp temperature inversion");
    require(std::abs(r.droplet_diameter_m-std::cbrt(6*m/(800*std::acos(-1.))))<1e-19,
            "mass-density diameter inversion");
  }
  auto cubic = a;
  cubic.pack.cp_j_per_kg_k = {TemperatureCorrelationKind::polynomial_cubic,300,{2000,2,.01,.0001}};
  const double theta = 40.;
  const double h = -1e6+2000*theta+theta*theta+.01/3*theta*theta*theta+
                   .0001/4*theta*theta*theta*theta;
  require(std::abs(remap_liquid_mass_enthalpy(cubic,m,h).temperature_k-340)<1e-11,
          "integrated cubic cp inverse");
  auto invalid = a;
  // cp positive at both endpoints and negative inside: reject multi-root h.
  invalid.pack.cp_j_per_kg_k = {TemperatureCorrelationKind::polynomial_cubic,400,{-1,0,1,0}};
  require(!remap_liquid_mass_enthalpy(invalid,m,-1e6).available,"interior cp minimum");
  invalid.pack.cp_j_per_kg_k = {TemperatureCorrelationKind::polynomial_cubic,400,{-1,0,1,.001}};
  require(!remap_liquid_mass_enthalpy(invalid,m,-1e6).available,"cubic derivative roots");
  for (double bad : {0.,-1.,std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()}) {
    const auto r = remap_liquid_mass_enthalpy(a,bad,-1e6);
    require(!r.available && r.temperature_k==0 && r.droplet_diameter_m==0,
            "invalid mass canonical failure");
  }
  for (double bad : {-1200001.,-239999.,std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()})
    require(!remap_liquid_mass_enthalpy(a,m,bad).available,"enthalpy bracket admission");
  a.pack.cp_j_per_kg_k.kind = TemperatureCorrelationKind::kerosene_cp_v1;
  a.pack.cp_j_per_kg_k.c = {};
  a.pack.cp_j_per_kg_k.reference_temperature_k = 298.15;
  a.pack.density_kg_per_m3.kind = TemperatureCorrelationKind::kerosene_density_v1;
  a.pack.density_kg_per_m3.c = {};
  a.pack.density_kg_per_m3.reference_temperature_k = 298.15;
  a.reference_temperature_k = 298.15;
  a.reference_liquid_enthalpy_j_per_kg = -1584371.4600473193;
  // Source NASA7 at its own R=8314.3 and the literal legacy density.
  const double coeff[]{2.0869217,.13314965,-8.1157452e-5,2.9409286e-8,-6.5195213e-12,-31310.966};
  for(double t : {351.2166748046875,400.,416.9096984863281}) {
    double vapor = coeff[5];
    for(unsigned i=0;i<5;++i) vapor += coeff[i]*std::pow(t,i+1)/(i+1);
    const double source_h = 8314.3/167.31771*vapor-
      250183*std::pow((684.26-t)/(684.26-483.15),.38);
    const auto r = remap_liquid_mass_enthalpy(a,m,source_h);
    require(r.available && r.temperature_k>t+.1 && r.temperature_k<t+.67,
            "approved source enthalpy drives target temperature shift");
    require(std::abs(r.enthalpy_residual_j_per_kg)<1e-8,"kerosene enthalpy conservation");
    const double oldrho = 1037.096-.7233865*t-9.255437+3/(733-t);
    const double old_d = std::cbrt(6*m/(std::acos(-1.)*oldrho));
    require(r.droplet_diameter_m>old_d*1.006 && r.droplet_diameter_m<old_d*1.010,
            "corrected density conserves source mass through diameter");
  }
  std::cout << "liquid remap analytic, nonmonotone, domain and actual-source checks passed\n";
}

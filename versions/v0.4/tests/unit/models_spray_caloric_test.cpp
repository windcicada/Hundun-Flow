// SPDX-License-Identifier: Apache-2.0
#include "models_chemistry_adapter_detail.hpp"
#include "models_spray_properties_detail.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <unistd.h>
using namespace hundun::v04;
using namespace spray::detail;
namespace {
bool close(double a, double b, double tolerance = 1e-11) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a-b) <= tolerance*std::max(1., std::abs(b));
}
std::string replace(std::string content, const char *key, const char *line) {
  const auto begin = content.find(std::string("\n")+key+' ')+1;
  content.replace(begin, content.find('\n', begin)-begin, line);
  return content;
}
// A thermodynamically consistent vapor shares h_v = h_l + L and
// cp_v = cp_l + dL/dT. This exercises the generic film's derivative audit.
class PhaseGas final : public portable::GasQueryProvider {
public:
  explicit PhaseGas(const LiquidAsset &asset) : asset_(asset) {}
  const portable::GasIdentity &gas_identity() const noexcept override {
    return asset_.gas_identity;
  }
  double latent(double t) const noexcept {
    return 250183*std::pow((684.26-t)/(684.26-483.15),.38);
  }
  double vapor_h(double t) const noexcept {
    return evaluate_liquid_enthalpy(asset_,t).liquid_enthalpy_j_per_kg+latent(t);
  }
  double mixture_h(double t, double y) const noexcept {
    return y*vapor_h(t)+(1-y)*1100*t;
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    double t = q.temperature_k;
    const double y = q.mass_fractions[0];
    if (q.coordinates == portable::GasStateCoordinates::pressure_enthalpy) {
      double lower=200, upper=680;
      for (unsigned i=0; i<64; ++i) {
        t=.5*(lower+upper);
        if (mixture_h(t,y)<q.enthalpy_j_per_kg) lower=t; else upper=t;
      }
      t=.5*(lower+upper);
    }
    const double cp = evaluate_liquid_enthalpy(asset_,t).cp_j_per_kg_k -
                       .38*latent(t)/(684.26-t);
    out.sample = {q.revision,q.composition_fingerprint,q.pressure_pa,t,
                  q.pressure_pa*28/(8314.46261815324*t),mixture_h(t,y),
                  y*cp+(1-y)*1100,2e-5,.05};
    out.species_enthalpies_j_per_kg[0]=vapor_h(t);
    out.species_enthalpies_j_per_kg[1]=1100*t;
    for (unsigned i=0; i<2; ++i) {
      out.diffusivities_m2_per_s[i]=2e-5;
      out.net_mass_rates_kg_per_m3_s[i]=0;
    }
    return portable::Status::success;
  }
private:
  const LiquidAsset &asset_;
};
}
int main(int argc, char **argv) {
  if (argc != 2) return 2;
  std::ifstream input(argv[1]);
  const std::string original((std::istreambuf_iterator<char>(input)), {});
  if (!input || original.empty()) return 2;
  char temporary[] = "/tmp/hf-cp-XXXXXX";
  const int fd = mkstemp(temporary);
  if (fd < 0) return 2;
  close(fd);
  struct Cleanup {const char *path; ~Cleanup(){std::filesystem::remove(path);}} cleanup{temporary};
  chemistry::detail::AnalyticIsomerBackend gas;
  auto load = [&](const std::string &content) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char c : content) {hash ^= c; hash *= UINT64_C(1099511628211);}
    {std::ofstream output(temporary); output << content;}
    return load_liquid_asset(temporary, hash, gas.gas_identity());
  };
  auto content = replace(original, "temperature_range", "temperature_range 200 680");
  content = replace(content, "density", "density kerosene_density_v1 298.15 0 0 0 0");
  content = replace(content, "cp", "cp kerosene_cp_v1 298.15 0 0 0 0");
  content = replace(content, "latent", "latent kerosene_latent_v1 298.15 0 0 0 0");
  content = replace(content, "saturation", "saturation kerosene_v1");
  auto loaded = load(content);
  if (!loaded.available) {
    std::cerr << "native liquid asset cannot carry THICK_EX caloric and saturation relations\n";
    return 1;
  }
  auto &asset = loaded.asset;
  const LiquidPropertyService service(&asset.pack, 1);
  // Independently evaluated with 70-digit decimal arithmetic from the
  // complete reference phase relation; h(298.15 K) = -100000 J/kg.
  const double reference[][5]{
      {200,1532.6308602108784,349367.79588083428,.013992118984083467,-277311.79425463732},
      {298.15,2067.5118052049447,320555.47709072102,186.97568861554933,-100000},
      {350,2319.1733337985552,303462.53276215075,2434.6433993167907,13815.197725137023},
      {477.95,2856.2052307690542,252621.73230145016,100000.00392026658,346097.44480411027},
      {600,3318.4414027092334,179758.38420133595,857130.91920754884,722319.53387475398},
      {680,6816.2628715143665,57826.732144478643,2242035.8377051679,1025222.1895430594}};
  bool passed = true;
  for (const auto &row : reference) {
    const auto p = service.evaluate({asset.content_fingerprint, row[0]});
    const auto h = evaluate_liquid_enthalpy(asset, row[0]);
    passed &= p.succeeded() && h.available && close(p.properties.cp_j_per_kg_k,row[1]) &&
        close(p.properties.latent_heat_j_per_kg,row[2]) &&
        close(p.properties.saturation_pressure_pa,row[3]) &&
        close(h.liquid_enthalpy_j_per_kg,row[4]) && close(h.cp_j_per_kg_k,row[1]);
  }
  for (double t : {250.,350.,477.95,600.,679.}) {
    const auto h0 = evaluate_liquid_enthalpy(asset,t);
    const auto hm = evaluate_liquid_enthalpy(asset,t-1e-3);
    const auto hp = evaluate_liquid_enthalpy(asset,t+1e-3);
    passed &= close((hp.liquid_enthalpy_j_per_kg-hm.liquid_enthalpy_j_per_kg)/.002,
                     h0.cp_j_per_kg_k,2e-8);
    auto shifted = asset;
    shifted.reference_temperature_k = t;
    shifted.reference_liquid_enthalpy_j_per_kg = h0.liquid_enthalpy_j_per_kg;
    for (const auto &row : reference)
      passed &= close(evaluate_liquid_enthalpy(shifted,row[0]).liquid_enthalpy_j_per_kg,row[4]);
  }
  auto zero = asset;
  zero.reference_liquid_enthalpy_j_per_kg = 0;
  const auto tiny = evaluate_liquid_enthalpy(zero,std::nextafter(298.15,400.));
  passed &= tiny.available && tiny.liquid_enthalpy_j_per_kg > 0 &&
      close(tiny.liquid_enthalpy_j_per_kg/(std::nextafter(298.15,400.)-298.15),reference[1][1]);
  passed &= !load(replace(content,"cp","cp kerosene_latent_v1 298.15 0 0 0 0")).available;
  passed &= !load(replace(content,"latent","latent kerosene_cp_v1 298.15 0 0 0 0")).available;
  passed &= !load(replace(content,"cp","cp kerosene_cp_v1 298.15 1 0 0 0")).available;
  passed &= !load(replace(content,"temperature_range","temperature_range 200 684.26")).available;
  passed &= !load(replace(content,"temperature_range","temperature_range 43 680")).available;
  passed &= !load(replace(content,"saturation","saturation kerosene_v1 1")).available;
  auto wrong = asset.pack;
  wrong.surface_tension_n_per_m = wrong.cp_j_per_kg_k;
  passed &= !LiquidPropertyService(&wrong,1).evaluate({asset.content_fingerprint,350}).succeeded();
  passed &= load(original).available;
  PhaseGas phase_gas(asset);
  FilmQueryWorkspace film(2);
  FilmQueryInput request;
  request.expected_revision={7,2,1};
  double y[]{.001,.999};
  request.surface_temperature_k=350;
  request.far_gas={request.expected_revision,asset.gas_identity.composition_fingerprint,
                  portable::GasStateCoordinates::pressure_enthalpy,101325,
                  phase_gas.mixture_h(400,y[0]),0,y,2};
  const auto sampled=film.query(asset,phase_gas,request);
  if (!sampled.available)
    std::cerr << "kerosene latent derivative fails shared film thermodynamic audit\n";
  passed &= sampled.available && std::abs(sampled.latent_heat_consistency_residual_j_per_kg)<1e-8;
  std::cout << "kerosene_caloric_asset_reference_domain passed=" << passed << '\n';
  return passed ? 0 : 1;
}

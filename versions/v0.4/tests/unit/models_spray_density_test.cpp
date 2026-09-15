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
using namespace hundun::v04::spray::detail;

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  // The new versioned density relation must be available through the same
  // provider used by parcel mass/diameter and transfer queries.
  constexpr auto density_kind = TemperatureCorrelationKind::kerosene_density_v1;
  LiquidPropertyPack pack;
  pack.material_fingerprint = 91;
  pack.minimum_temperature_k = 280;
  pack.maximum_temperature_k = 680;
  pack.density_kg_per_m3 = {density_kind, 298.15, {}};
  pack.cp_j_per_kg_k.c[0] = 2200;
  pack.latent_heat_j_per_kg.c[0] = 250000;
  pack.surface_tension_n_per_m.c[0] = .025;
  pack.viscosity_pa_s.c[0] = .001;
  pack.saturation_pressure.antoine_a = 4;
  pack.saturation_pressure.pressure_scale_pa = 1;
  const LiquidPropertyService service(&pack, 1);
  // Decimal evaluations of Rachner (1998), equations 3.10.3/3.10.4, p.48.
  const double reference[][2]{{298.15,800.1341089769375},
      {350,759.745093146214}, {400,719.9472948948949}, {500,635.6798444206008}};
  bool passed = true;
  for (const auto& row : reference) {
    const auto r = service.evaluate({91, row[0]});
    passed &= r.succeeded() && std::abs(r.properties.density_kg_per_m3-row[1]) < 1e-10;
    const auto phase = evaluate_kerosene_phase(row[0],101325.,298.15);
    passed &= phase.available && phase.liquid_density_kg_per_m3 == r.properties.density_kg_per_m3;
  }
  if (!passed) {
    std::cerr << "kerosene density provider fails published coefficient samples\n";
    return 1;
  }
  double previous = std::numeric_limits<double>::infinity();
  for (double t = 280; t <= 680; t += 2) {
    const auto r = service.evaluate({91,t});
    const double rho = r.properties.density_kg_per_m3;
    passed &= r.succeeded() && rho > 0 && rho < previous;
    previous = rho;
  }
  auto bad = pack;
  bad.maximum_temperature_k = 800;
  const LiquidPropertyService out_of_domain(&bad,1);
  for (double t : {684.26,733.,800.,std::numeric_limits<double>::quiet_NaN()})
    passed &= !out_of_domain.evaluate({91,t}).succeeded();
  bad = pack; bad.density_kg_per_m3.c[0] = 1;
  passed &= !LiquidPropertyService(&bad,1).evaluate({91,350}).succeeded();
  bad = pack; bad.cp_j_per_kg_k = pack.density_kg_per_m3;
  passed &= !LiquidPropertyService(&bad,1).evaluate({91,350}).succeeded();

  std::ifstream input(argv[1]);
  const std::string original((std::istreambuf_iterator<char>(input)), {});
  if (!input || original.empty()) return 2;
  char temporary[] = "/tmp/hf-rho-XXXXXX";
  const int descriptor = mkstemp(temporary);
  if (descriptor < 0) return 2;
  close(descriptor);
  struct Cleanup {const char* path; ~Cleanup(){std::filesystem::remove(path);}} cleanup{temporary};
  chemistry::detail::AnalyticIsomerBackend gas;
  const auto load = [&](std::string content) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char c : content) {hash ^= c; hash *= UINT64_C(1099511628211);}
    {std::ofstream file(temporary); file << content;}
    return load_liquid_asset(temporary, hash, gas.gas_identity());
  };
  const auto replace = [](std::string content,const char* name,const char* line) {
    const auto begin = content.find(std::string("\n") + name + ' ') + 1;
    content.replace(begin,content.find('\n',begin)-begin,line);
    return content;
  };
  const auto content = replace(original,"density","density kerosene_density_v1 298.15 0 0 0 0");
  const auto loaded = load(content);
  passed &= loaded.available && loaded.asset.pack.material_fingerprint == loaded.asset.content_fingerprint;
  if (loaded.available) {
    const auto r = LiquidPropertyService(&loaded.asset.pack,1).evaluate({loaded.asset.content_fingerprint,350});
    passed &= r.succeeded() && std::abs(r.properties.density_kg_per_m3-759.745093146214)<1e-10;
  }
  passed &= !load(replace(content,"cp","cp kerosene_density_v1 298.15 0 0 0 0")).available;
  passed &= !load(replace(content,"density","density kerosene_density_v1 298.15 1 0 0 0")).available;
  passed &= !load(replace(content,"temperature_range","temperature_range 280 700")).available;
  const auto old = load(original);
  passed &= old.available && old.asset.content_fingerprint != loaded.asset.content_fingerprint;
  std::cout << "kerosene_density published_values_domain_asset_identity passed=" << passed << '\n';
  return passed ? 0 : 1;
}

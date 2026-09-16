// SPDX-License-Identifier: Apache-2.0
// Offline SI mass/enthalpy table adapter. The JSONL orchestrator retains
// source identity and history, validates the whole result, then publishes it.
#include "hundun/v04_cantera.hpp"
#include "models_spray_remap_detail.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace hundun::v04;
int main(int argc, char** argv) {
  try {
    if (argc != 7) throw std::runtime_error("usage: v04_spray_remap GAS SHA PHASE LIQUID FNV SPECIES_CSV < mass_h.dat");
    chemistry::CanteraBackendConfig config;
    config.mechanism = {argv[1],argv[2],argv[3]};
    config.chemistry = {1e-10,1e-16,20000};
    config.continuous_enthalpy = true;
    std::istringstream names(argv[6]);
    std::string name;
    while (std::getline(names,name,',')) config.species_names.push_back(name);
    auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config);
    chemistry::CanteraWorkspacePool pool(runtime,1);
    auto backend = chemistry::make_cantera_backend(config,pool);
    std::size_t used = 0;
    const std::string expected(argv[5]);
    const auto fnv = std::stoull(expected,&used);
    if (used != expected.size()) throw std::runtime_error("invalid liquid identity");
    const auto asset = spray::detail::load_liquid_asset(argv[4],fnv,backend->gas_identity());
    if (!asset.available) throw std::runtime_error("liquid/gas asset admission failed");
    std::cout << std::setprecision(17);
    std::string line;
    std::size_t row = 0;
    while (std::getline(std::cin,line)) {
      ++row;
      std::istringstream in(line);
      double mass{}, h{};
      std::string extra;
      if (!(in >> mass >> h) || (in >> extra))
        throw std::runtime_error("invalid numeric row "+std::to_string(row));
      const auto r = spray::detail::remap_liquid_mass_enthalpy(asset.asset,mass,h);
      if (!r.available) throw std::runtime_error("liquid remap failed at row "+std::to_string(row));
      std::cout << r.temperature_k << ' ' << r.droplet_diameter_m << ' '
        << r.density_kg_per_m3 << ' ' << r.specific_enthalpy_j_per_kg << ' '
        << r.enthalpy_residual_j_per_kg << ' ' << r.geometry_relative_mass_error
        << ' ' << asset.asset.pack.material_fingerprint << '\n';
    }
    if (!std::cin.eof() || row==0 || !std::cout)
      throw std::runtime_error("empty or failed table stream");
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_combustion.hpp"

#include <cmath>
#include <iostream>

using namespace hundun::v04::combustion;
int main() {
  bool passed = true;
  auto expect = [&](bool value, const char* message) {
    if (!value) std::cerr << message << '\n';
    passed &= value;
  };
  ChemistryIdentity identity;
  identity.element_count = 1U;
  identity.species = {{1.0, {1U}, 0.0}, {1.0, {1U}, 0.0}};
  identity.fingerprint = chemistry_identity_fingerprint(identity);
  ThermochemicalState mean{101325.0, 1.0, 200.0, {0.5, 0.5}};
  const ThermochemicalState lean{101325.0, 1.0, 100.0, {0.2, 0.8}};
  const ThermochemicalState rich{101325.0, 1.0, 300.0, {0.8, 0.2}};
  EsfCommonSource source;
  source.composition_fingerprint = identity.fingerprint;
  source.density_delta_kg_per_m3 = 0.25;
  source.species_density_delta_kg_per_m3 = {0.25, 0.0};
  source.thermochemical_enthalpy_density_delta_j_per_m3 = 50.0;
  for (const std::size_t count : {2U, 4U}) {
    std::vector<ThermochemicalState> fields;
    for (std::size_t i = 0; i < count; ++i) fields.push_back(i % 2U ? rich : lean);
    const auto result = apply_esf_common_source(mean, fields, identity, source);
    expect(result.succeeded() && result.fields.size() == count &&
           result.model_id == "esf_common_conservative_source_v1",
           "one common integrated source creates all field candidates");
    if (!result.succeeded()) continue;
    expect(result.mean_state.density_kg_per_m3 == 1.25 &&
           std::abs(result.mean_state.mass_fractions[0] - 0.6) < 1e-14 &&
           result.mean_state.total_thermochemical_enthalpy_j_per_kg == 200.0,
           "mean inventory receives source once, independent of N");
    expect(std::abs(result.fields[0].mass_fractions[0] - 0.36) < 1e-14 &&
           std::abs(result.fields[1].mass_fractions[0] - 0.84) < 1e-14 &&
           result.fields[0].total_thermochemical_enthalpy_j_per_kg == 120.0 &&
           result.fields[1].total_thermochemical_enthalpy_j_per_kg == 280.0,
           "each field receives the same physical source, not source divided by N");
    expect(fields[0].mass_fractions == lean.mass_fractions &&
           mean.mass_fractions[0] == 0.5,
           "mean and all input fields remain unchanged");
    double average = 0.0, variance = 0.0;
    for (const auto& field : result.fields) {
      average += field.mass_fractions[0] / static_cast<double>(count);
      const double deviation = field.mass_fractions[0] - 0.6;
      variance += deviation * deviation / static_cast<double>(count);
    }
    expect(std::abs(average - 0.6) < 1e-14 && std::abs(variance - 0.0576) < 1e-14,
           "mean and scalar variance follow the conservative common-source map");
  }
  std::vector<ThermochemicalState> fields{lean, rich};
  source.species_density_delta_kg_per_m3 = {-0.3, 0.55};
  auto failed = apply_esf_common_source(mean, fields, identity, source);
  expect(!failed.succeeded() && !failed.available && failed.fields.empty() &&
         failed.mean_state.mass_fractions.empty(),
         "one inadmissible field rejects the whole candidate, without clipping");
  source.species_density_delta_kg_per_m3 = {0.5, 0.0};
  failed = apply_esf_common_source(mean, fields, identity, source);
  expect(!failed.succeeded() && failed.fields.empty(),
         "species source sum must match total mass source");
  source.species_density_delta_kg_per_m3 = {0.25, 0.0};
  ++source.composition_fingerprint;
  failed = apply_esf_common_source(mean, fields, identity, source);
  expect(!failed.succeeded() && failed.fields.empty(), "source identity mismatch is explicit");
  --source.composition_fingerprint;
  mean.mass_fractions = {0.6, 0.4};
  failed = apply_esf_common_source(mean, fields, identity, source);
  expect(!failed.succeeded() && failed.fields.empty(),
         "inconsistent initial mean is not silently replaced");
  return passed ? 0 : 1;
}

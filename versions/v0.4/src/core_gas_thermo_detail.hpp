// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_physics.hpp"
#include "hundun/v04_portable.hpp"
#include "physics_input_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::detail {
// Cold compatibility check of model definitions, never a fit to sampled
// states. Allow only floating-point serialization roundoff in coefficients;
// an absolute floor would incorrectly admit altered tiny high-order terms.
inline bool native_nasa7_compatible(const ThermophysicalSpec& input,
    const portable::GasIdentity& gas, portable::IdealGasNasa7View provider) {
  if (!provider.species || provider.species_count!=input.species.size() ||
      provider.species_count!=gas.species_names.size() ||
      provider.species_count!=gas.molecular_weights_kg_per_kmol.size() ||
      !gas.composition_fingerprint ||
      provider.composition_fingerprint!=gas.composition_fingerprint ||
      provider.gas_constant_j_per_kmol_k!=kUniversalGasConstant ||
      !(provider.minimum_temperature_k>0) ||
      !std::isfinite(provider.minimum_temperature_k) ||
      !std::isfinite(provider.maximum_temperature_k) ||
      provider.minimum_temperature_k>input.minimum_temperature ||
      provider.maximum_temperature_k<input.maximum_temperature) return false;
  auto native=input;
  if(!canonicalize_thermophysical_spec(native))return false;
  const auto same=[](double a,double b) {
    return std::isfinite(a) && std::isfinite(b) &&
        (a==b || std::abs(a-b)<=64*std::numeric_limits<double>::epsilon()*
                                      std::max(std::abs(a),std::abs(b)));
  };
  for(std::size_t i=0;i<provider.species_count;++i) {
    for(std::size_t j=0;j<i;++j)
      if(gas.species_names[i]==gas.species_names[j])return false;
    const auto found=std::find_if(native.species.begin(),native.species.end(),
        [&](const auto& s){return s.stable_name==gas.species_names[i];});
    if(found==native.species.end() || found->transport_law==TransportLaw::nasa_air)return false;
    const auto& supplied=provider.species[i];
    if(!same(found->molecular_weight,supplied.molecular_weight_kg_per_kmol) ||
        !same(gas.molecular_weights_kg_per_kmol[i],supplied.molecular_weight_kg_per_kmol) ||
        found->temperature_switch!=supplied.temperature_switch_k)return false;
    for(std::size_t c=0;c<7;++c)
      if(!same(found->nasa7_low[c],supplied.low[c]) ||
          !same(found->nasa7_high[c],supplied.high[c]))return false;
  }
  return true;
}
} // namespace hundun::v04::detail

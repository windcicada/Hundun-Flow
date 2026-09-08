// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_case.hpp"
#include <cmath>
#include <limits>

namespace hundun::v04::detail {
inline bool valid_spray_spec(const SpraySpec &s) {
  const auto finite = [](Real3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  const auto positive = [](double x) { return std::isfinite(x) && x > 0; };
  if (s.liquid_file.empty() || s.liquid_file.has_parent_path() ||
      s.liquid_file.extension() != ".asset" ||
      s.liquid_file.native().size() > 255 ||
      s.liquid_file.native().find('\0') != std::string::npos ||
      !s.liquid_fingerprint || !s.maximum_local_parcels ||
      s.maximum_local_parcels >
          std::uint32_t(std::numeric_limits<int>::max() / 18) ||
      !s.maximum_local_segments ||
      s.maximum_local_segments >
          std::uint32_t(std::numeric_limits<int>::max() / 32) ||
      !positive(s.maximum_substep_s) || !positive(s.minimum_substep_s) ||
      s.minimum_substep_s > s.maximum_substep_s ||
      !positive(s.relative_tolerance) || s.relative_tolerance >= 1 ||
      s.injectors.empty() || s.injectors.size() > 64)
    return false;
  for (std::size_t i = 0; i < s.injectors.size(); ++i) {
    const auto &v = s.injectors[i];
    const double axis_norm = std::hypot(v.axis.x, v.axis.y, v.axis.z);
    if (!v.id || !finite(v.origin_m) || !finite(v.axis) ||
        !positive(axis_norm) || !std::isfinite(v.cone_half_angle_rad) ||
        v.cone_half_angle_rad < 0 ||
        v.cone_half_angle_rad >= 1.5707963267948966 ||
        !std::isfinite(v.speed_m_per_s) || v.speed_m_per_s < 0 ||
        !std::isfinite(v.mass_flow_rate_kg_per_s) ||
        v.mass_flow_rate_kg_per_s < 0 ||
        !positive(v.represented_mass_per_parcel_kg) ||
        !positive(v.droplet_diameter_m) || !positive(v.temperature_k))
      return false;
    for (std::size_t j = 0; j < i; ++j)
      if (s.injectors[j].id == v.id)
        return false;
  }
  return true;
}

template <class Writer> bool write_spray(Writer &w, const SpraySpec &s) {
  if (!w.text(s.liquid_file.generic_string()))
    return false;
  w.u64(s.liquid_fingerprint);
  w.u64(s.seed);
  w.u32(s.maximum_local_parcels);
  w.u32(s.maximum_local_segments);
  w.real(s.maximum_substep_s);
  w.real(s.minimum_substep_s);
  w.real(s.relative_tolerance);
  w.byte(s.tab_breakup);
  w.u32(static_cast<std::uint32_t>(s.injectors.size()));
  for (const auto &v : s.injectors) {
    w.u64(v.id);
    w.real3(v.origin_m);
    w.real3(v.axis);
    w.real(v.cone_half_angle_rad);
    w.real(v.speed_m_per_s);
    w.real(v.mass_flow_rate_kg_per_s);
    w.real(v.represented_mass_per_parcel_kg);
    w.real(v.droplet_diameter_m);
    w.real(v.temperature_k);
  }
  return true;
}

template <class Reader> bool read_spray(Reader &r, SpraySpec &s) {
  std::string file;
  std::uint8_t tab{};
  std::uint32_t count{};
  if (!r.text(file) || !r.u64(s.liquid_fingerprint) || !r.u64(s.seed) ||
      !r.u32(s.maximum_local_parcels) || !r.u32(s.maximum_local_segments) ||
      !r.real(s.maximum_substep_s) || !r.real(s.minimum_substep_s) ||
      !r.real(s.relative_tolerance) || !r.byte(tab) || tab > 1 ||
      !r.u32(count) || count == 0 || count > 64)
    return false;
  s.liquid_file = std::move(file);
  s.tab_breakup = tab != 0;
  s.injectors.resize(count);
  for (auto &v : s.injectors)
    if (!r.u64(v.id) || !r.real3(v.origin_m) || !r.real3(v.axis) ||
        !r.real(v.cone_half_angle_rad) || !r.real(v.speed_m_per_s) ||
        !r.real(v.mass_flow_rate_kg_per_s) ||
        !r.real(v.represented_mass_per_parcel_kg) ||
        !r.real(v.droplet_diameter_m) || !r.real(v.temperature_k))
      return false;
  return valid_spray_spec(s);
}

template <class Hash> void hash_spray(Hash &h, const SpraySpec &s) {
  h.text("spray-case-v1");
  h.text(s.liquid_file.generic_string());
  h.integer(s.liquid_fingerprint);
  h.integer(s.seed);
  h.integer(s.maximum_local_parcels);
  h.integer(s.maximum_local_segments);
  h.real(s.maximum_substep_s);
  h.real(s.minimum_substep_s);
  h.real(s.relative_tolerance);
  h.integer(static_cast<std::uint8_t>(s.tab_breakup));
  h.integer(s.injectors.size());
  for (const auto &v : s.injectors) {
    h.integer(v.id);
    for (double x :
         {v.origin_m.x, v.origin_m.y, v.origin_m.z, v.axis.x, v.axis.y,
          v.axis.z, v.cone_half_angle_rad, v.speed_m_per_s,
          v.mass_flow_rate_kg_per_s, v.represented_mass_per_parcel_kg,
          v.droplet_diameter_m, v.temperature_k})
      h.real(x);
  }
}
} // namespace hundun::v04::detail

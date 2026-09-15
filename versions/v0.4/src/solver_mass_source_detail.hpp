// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"
namespace hundun::v04::detail {
// Equivalent conservative mass of the frozen-density advective scalar row.
// Callers validate views/a0 before entry and validate the resulting mass.
inline double frozen_density_carrier_mass(ConstFaceFluxView flux, Int3 c,
                                          double accepted_mass,
                                          double a0) noexcept {
  const long double net=(static_cast<long double>(flux.x.unchecked({c.x+1,c.y,c.z}))-flux.x.unchecked(c))+
      (static_cast<long double>(flux.y.unchecked({c.x,c.y+1,c.z}))-flux.y.unchecked(c))+
      (static_cast<long double>(flux.z.unchecked({c.x,c.y,c.z+1}))-flux.z.unchecked(c));
  return static_cast<double>(accepted_mass-net/a0);
}
inline bool valid_mass_source(ConservativeMassSourceView source,
                              PlanFingerprint identity, RevisionToken time,
                              Int3 cells) noexcept {
  if (identity == 0U) {
    return source.identity == 0U && source.time == 0U &&
           source.rate.base == nullptr && source.rate.revision == 0U &&
           source.rate.storage_identity == 0U &&
           source.rate.revision_domain == 0U;
  }
  return source.identity == identity && time != 0U && source.time == time &&
         valid_cell_view(source.rate, cells, 0U, 1U, 0U);
}
inline double mass_source_rate(ConservativeMassSourceView source,
                               Int3 cell) noexcept {
  return source.identity == 0U ? 0.0 : source.rate.unchecked(cell, 0U);
}
inline bool same_mass_source(ConservativeMassSourceView a,
                             ConservativeMassSourceView b) noexcept {
  return a.identity == b.identity && a.time == b.time &&
         a.rate.base == b.rate.base && a.rate.field == b.rate.field &&
         a.rate.replica == b.rate.replica &&
         a.rate.revision == b.rate.revision &&
         a.rate.storage_identity == b.rate.storage_identity &&
         a.rate.revision_domain == b.rate.revision_domain;
}
} // namespace hundun::v04::detail

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_spray.hpp"
#include <cmath>
#include <cstring>

namespace hundun::v04::spray::detail {
// Version zero carries the canonical empty state. Version one stores SI
// unit-mass dissipation and independently accumulated exposure/rate clocks.
struct PersistentSgsBreakupState {
  std::uint64_t version{};
  SgsBreakupHistory history{};
};
inline constexpr std::size_t kSgsHistoryLanes = 6U;

inline bool valid_sgs_state(const PersistentSgsBreakupState& state) noexcept {
  const auto& h = state.history;
  if (state.version > 1U || h.poisson_multiplier > 7U) return false;
  for (double value : {h.mean_dissipation_m2_per_s3, h.dissipation_age_s,
                      h.mean_rate_per_s, h.rate_age_s}) {
    if (!std::isfinite(value) || value < 0.0 || (state.version == 0U && value != 0.0))
      return false;
  }
  return state.version == 1U || h.poisson_multiplier == 0U;
}
inline void encode_sgs_state(const PersistentSgsBreakupState& state,
                             std::uint64_t* wire) noexcept {
  wire[0] = state.version;
  const auto& h = state.history;
  const double values[]{h.mean_dissipation_m2_per_s3, h.dissipation_age_s,
                        h.mean_rate_per_s, h.rate_age_s};
  for (std::size_t i=0;i<4U;++i) std::memcpy(wire+1U+i,values+i,8U);
  wire[5] = h.poisson_multiplier;
}
inline PersistentSgsBreakupState decode_sgs_state(const std::uint64_t* wire) noexcept {
  PersistentSgsBreakupState state;
  state.version = wire[0];
  double values[4]{};
  for (std::size_t i=0;i<4U;++i) std::memcpy(values+i,wire+1U+i,8U);
  state.history = {values[0],values[1],values[2],values[3],
      static_cast<std::uint8_t>(wire[5] <= 7U ? wire[5] : 8U)};
  return state;
}
} // namespace hundun::v04::spray::detail

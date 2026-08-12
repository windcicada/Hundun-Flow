// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/les_wale.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace hundun::flow::detail {

class ImmersedWaleAttemptAuthority final {
public:
  ImmersedWaleAttemptAuthority(
      les::WaleAttemptCoefficients coefficients,
      std::vector<double> effective_dynamic_viscosity_by_active_cell,
      std::vector<double> effective_dynamic_viscosity_by_face,
      std::uint64_t wall_effective_viscosity_fingerprint)
      : coefficients_(std::move(coefficients)),
        effective_dynamic_viscosity_by_active_cell_(
            std::move(effective_dynamic_viscosity_by_active_cell)),
        effective_dynamic_viscosity_by_face_(
            std::move(effective_dynamic_viscosity_by_face)),
        summary_(coefficients_.summary()),
        wall_effective_viscosity_fingerprint_(
            wall_effective_viscosity_fingerprint) {}

  ~ImmersedWaleAttemptAuthority() noexcept = default;
  ImmersedWaleAttemptAuthority(ImmersedWaleAttemptAuthority &&) noexcept =
      default;
  ImmersedWaleAttemptAuthority &
  operator=(ImmersedWaleAttemptAuthority &&) = delete;
  ImmersedWaleAttemptAuthority(const ImmersedWaleAttemptAuthority &) = delete;
  ImmersedWaleAttemptAuthority &
  operator=(const ImmersedWaleAttemptAuthority &) = delete;

  const std::vector<double> &
  effective_dynamic_viscosity_by_active_cell() const noexcept {
    return effective_dynamic_viscosity_by_active_cell_;
  }

  const std::vector<double> &
  effective_dynamic_viscosity_by_face() const noexcept {
    return effective_dynamic_viscosity_by_face_;
  }

  const les::WaleSummary &summary() const noexcept { return summary_; }

  std::uint64_t evaluation_count() const noexcept { return 1U; }

  std::uint64_t wall_effective_viscosity_fingerprint() const noexcept {
    return wall_effective_viscosity_fingerprint_;
  }

private:
  les::WaleAttemptCoefficients coefficients_;
  std::vector<double> effective_dynamic_viscosity_by_active_cell_;
  std::vector<double> effective_dynamic_viscosity_by_face_;
  les::WaleSummary summary_;
  std::uint64_t wall_effective_viscosity_fingerprint_{};
};

} // namespace hundun::flow::detail

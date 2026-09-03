// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/exec_execution.hpp"
#include "hundun/diag_stage3_performance.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_view.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::les {

enum class WaleTimeOrder : std::uint8_t { backward_euler = 1, bdf2 = 2 };

struct WaleControl final {
  double coefficient{};
  double turbulent_prandtl{};
  double turbulent_schmidt{};
};

struct WaleCoefficientIdentity final {
  std::uint64_t value{};
  friend bool operator==(WaleCoefficientIdentity left,
                         WaleCoefficientIdentity right) noexcept {
    return left.value == right.value;
  }
  friend bool operator!=(WaleCoefficientIdentity left,
                         WaleCoefficientIdentity right) noexcept {
    return !(left == right);
  }
};

struct WaleAttemptInput final {
  std::uint64_t step{};
  double attempted_dt_s{};
  WaleTimeOrder order{WaleTimeOrder::backward_euler};
  std::uint64_t committed_state_fingerprint{};
  std::uint64_t history_state_fingerprint{};
  std::uint64_t lagged_gradient_fingerprint{};
  std::uint64_t density_fingerprint{};
  runtime::FieldView<const double> lagged_velocity_gradient;
  runtime::FieldView<const double> rho_attempt;
};

struct WaleSummary final {
  WaleCoefficientIdentity identity{};
  double minimum_nu_t_m2_per_s{};
  double maximum_nu_t_m2_per_s{};
  double l2_nu_t_m2_per_s{};
  std::uint64_t exact_zero_count{};
  std::uint64_t owned_active_count{};
};

class WaleAttemptCoefficients final {
public:
  ~WaleAttemptCoefficients() noexcept;
  WaleAttemptCoefficients(WaleAttemptCoefficients &&) noexcept;
  WaleAttemptCoefficients &operator=(WaleAttemptCoefficients &&) = delete;
  WaleAttemptCoefficients(const WaleAttemptCoefficients &) = delete;
  WaleAttemptCoefficients &operator=(const WaleAttemptCoefficients &) = delete;

  WaleCoefficientIdentity identity() const noexcept;
  const WaleSummary &summary() const noexcept;
  std::size_t owned_active_count() const noexcept;
  std::size_t local_active_count() const noexcept;
  execution::VectorView<const double> nu_t_m2_per_s() const;
  execution::VectorView<const double> mu_sgs_pa_s() const;

private:
  struct Impl;
  explicit WaleAttemptCoefficients(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class WaleModel;
};

class WaleModel final {
public:
  static WaleModel create(
      WaleControl, const mesh::MeshTopology &, const mesh::MeshGeometry &,
      std::size_t owned_active_count,
      const std::vector<mesh::GlobalCellId> &ordered_local_active_global_cells,
      execution::ExecutionContext &);

  ~WaleModel() noexcept;
  WaleModel(WaleModel &&) noexcept;
  WaleModel &operator=(WaleModel &&) = delete;
  WaleModel(const WaleModel &) = delete;
  WaleModel &operator=(const WaleModel &) = delete;

  WaleAttemptCoefficients evaluate(const WaleAttemptInput &) const;
  WaleControl control() const noexcept;
  diagnostics::Stage3PerformanceCounters
  performance_counters() const noexcept;

private:
  struct Impl;
  explicit WaleModel(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
};

} // namespace hundun::les

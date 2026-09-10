// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "solver_cartesian_detail.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace hundun::v04::detail {
// Initial guesses for the same correction stage on the immediately preceding
// composition sweep. They carry no operator, residual or convergence authority.
// Krylov always evaluates the current true residual and its ordinary tolerance.
class PressureEnergyDirectionHistory {
public:
  void reserve(Int3 shape) {
    shape_ = shape;
    const auto count = std::size_t(shape.x) * shape.y * shape.z;
    for (auto &values : values_)
      values.resize(count);
    sweeps_.fill(0U);
    sweep_ = 0U;
    previous_residual_ = std::numeric_limits<double>::infinity();
    restarts_ = 0U;
  }
  void begin_coupling(unsigned sweep) noexcept {
    // A new proposal or dt retry starts at sweep one. Skipped/reversed sweep
    // sequences cannot borrow an older proposal's direction.
    if (sweep <= 1U || sweep != sweep_ + 1U) {
      sweeps_.fill(0U);
      previous_residual_ = std::numeric_limits<double>::infinity();
    }
    sweep_ = sweep;
  }
  void observe_coupling_residual(double residual) noexcept {
    // A warm direction may already satisfy the linear absolute tolerance
    // while retaining an error above the nonlinear composition tolerance.
    // Once progress stalls, the next sweep uses the ordinary solver seed.
    // This changes neither equation residuals nor their acceptance gates.
    if (!std::isfinite(residual) || residual < 0.0 ||
        residual >= 0.9 * previous_residual_) {
      sweeps_.fill(0U);
      ++restarts_;
    }
    previous_residual_ = residual;
  }
  std::uint64_t restarts() const noexcept { return restarts_; }
  Status restore(unsigned corrector, unsigned refinement,
                 FieldView out) const noexcept {
    const int slot = select(corrector, refinement);
    if (slot < 0 || sweep_ <= 1U || sweeps_[slot] != sweep_ - 1U)
      return {};
    if (!valid_cell_view(out, shape_, 0U, 1U))
      return {StatusCode::invalid_plan, 1585U};
    each([&](Int3 c, std::size_t i) {
      out.unchecked(c, 0U) = values_[slot][i];
    });
    return {};
  }
  Status remember(unsigned corrector, unsigned refinement,
                  ConstFieldView in) noexcept {
    const int slot = select(corrector, refinement);
    if (slot < 0)
      return {};
    if (sweep_ == 0U || !valid_cell_view(in, shape_, 0U, 1U, 0U))
      return {StatusCode::invalid_plan, 1585U};
    bool finite = true;
    each([&](Int3 c, std::size_t) {
      finite &= std::isfinite(in.unchecked(c, 0U));
    });
    if (!finite)
      return {StatusCode::numerical_failure, 1585U};
    each(
        [&](Int3 c, std::size_t i) { values_[slot][i] = in.unchecked(c, 0U); });
    sweeps_[slot] = sweep_;
    return {};
  }
  std::uint64_t owned_payload_bytes() const noexcept {
    std::uint64_t bytes = 0U;
    for (const auto &values : values_)
      bytes += values.capacity() * sizeof(double);
    return bytes;
  }

private:
  static int select(unsigned corrector, unsigned refinement) noexcept {
    if (corrector == 1U && refinement == 0U)
      return 0;
    if (corrector == 2U && refinement <= 2U)
      return int(refinement) + 1;
    return -1;
  }
  template <class F> void each(F f) const noexcept {
    std::size_t i = 0U;
    for (int z = 0; z < shape_.z; ++z)
      for (int y = 0; y < shape_.y; ++y)
        for (int x = 0; x < shape_.x; ++x, ++i)
          f(Int3{x, y, z}, i);
  }
  Int3 shape_{};
  unsigned sweep_{};
  double previous_residual_{std::numeric_limits<double>::infinity()};
  std::uint64_t restarts_{};
  std::array<unsigned, 4> sweeps_{};
  std::array<std::vector<double>, 4> values_;
};
} // namespace hundun::v04::detail

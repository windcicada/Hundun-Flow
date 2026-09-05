// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

// Disjoint local elapsed intervals. No allocation, synchronization, or MPI;
// unwinding and early returns account for the active phase as well.
template <std::size_t N>
class LocalPhaseTimer {
 public:
  explicit LocalPhaseTimer(std::array<std::uint64_t, N>& totals) noexcept
      : totals_(totals), last_(Clock::now()) {}
  ~LocalPhaseTimer() noexcept { stop(); }
  LocalPhaseTimer(const LocalPhaseTimer&) = delete;
  LocalPhaseTimer& operator=(const LocalPhaseTimer&) = delete;
  void phase(std::size_t next) noexcept {
    if (next >= N || stopped_) return;
    account();
    current_ = next;
  }
  void stop() noexcept {
    if (!stopped_) account();
    stopped_ = true;
  }

 private:
  using Clock = std::chrono::steady_clock;
  void account() noexcept {
    const auto now = Clock::now();
    const auto elapsed = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - last_)
            .count());
    auto& total = totals_[current_];
    total = elapsed > UINT64_MAX - total ? UINT64_MAX : total + elapsed;
    last_ = now;
  }
  std::array<std::uint64_t, N>& totals_;
  Clock::time_point last_;
  std::size_t current_{};
  bool stopped_{};
};

}  // namespace hundun::v04::detail

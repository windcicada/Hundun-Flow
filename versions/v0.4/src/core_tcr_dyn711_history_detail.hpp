// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_io.hpp"
#include "models_tcr_dyn711_detail.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace hundun::v04::detail {
// Independent window identity: 624CF per-step rates cannot initialize dyn711
// cumulative rates. Each cell record carries the same accepted model clocks,
// allowing the generic Restart owner to redistribute complete records.
class Dyn711History {
public:
  using Rate = tcr::detail::Dyn711RateState;
  using Clock = tcr::detail::Dyn711Clock;
  void configure(PlanFingerprint model, std::size_t cells, std::size_t species,
                 double initial_cphi) {
    ns_ = species;
    width_ = static_cast<std::uint32_t>(40 + 40 * ns_);
    identity_ = (model ^ UINT64_C(0x5443523731310001) ^ ns_) * UINT64_C(1099511628211);
    if (!identity_) identity_ = 1;
    step_ = calls_ = 0;
    discard();
    accepted_.assign(cells * ns_, Rate{});
    trial_ = accepted_;
    cphi_.assign(cells, initial_cphi);
    trial_cphi_ = cphi_;
    staged_.resize(cells * ns_);
    cphi_staged_.resize(cells);
    inactive_.resize(cells);
    bytes_.resize(cells * width_);
    trial_bytes_.resize(bytes_.size());
    encode(accepted_, cphi_, 0, 0, bytes_);
  }
  std::uint64_t owned_bytes() const noexcept {
    return sizeof(*this) + sizeof(Rate) * (accepted_.capacity() + trial_.capacity()) +
        8 * (cphi_.capacity() + trial_cphi_.capacity()) + staged_.capacity() +
        cphi_staged_.capacity() + inactive_.capacity() + bytes_.capacity() + trial_bytes_.capacity();
  }
  const Rate &accepted(std::size_t cell, std::size_t species) const noexcept {
    return accepted_[cell * ns_ + species];
  }
  std::size_t cells() const noexcept { return cphi_.size(); }
  std::size_t species() const noexcept { return ns_; }
  double cphi(std::size_t cell) const noexcept { return cphi_[cell]; }
  Clock clock() const noexcept { return clock_at(calls_); }
  std::uint64_t statistics_calls() const noexcept { return calls_; }
  Status begin(std::uint64_t step) noexcept {
    discard();
    if (step != step_ || step == UINT64_MAX || calls_ == UINT64_MAX || ns_ == 0)
      return invalid();
    std::copy(accepted_.begin(), accepted_.end(), trial_.begin());
    std::copy(cphi_.begin(), cphi_.end(), trial_cphi_.begin());
    std::fill(staged_.begin(), staged_.end(), 0);
    std::fill(inactive_.begin(), inactive_.end(), 0);
    std::fill(cphi_staged_.begin(), cphi_staged_.end(),
        tcr::detail::dyn711_tick(clock()).update_cphi ? 0 : 1);
    candidate_step_ = step + 1;
    candidate_calls_ = calls_ + 1;
    preparing_ = true;
    return {};
  }
  Status stage_rate(std::size_t cell, std::size_t species, double dt,
      double pdf_rate, double psr_rate, double eta, double chemical_time,
      double flow_time, double weak_sum) noexcept {
    if (!preparing_ || cell >= cphi_.size() || species >= ns_ || inactive_[cell])
      return invalid();
    const auto i = cell * ns_ + species;
    if (!tcr::detail::dyn711_advance_rate(accepted_[i], clock(), dt, pdf_rate, psr_rate,
        eta, chemical_time, flow_time, weak_sum, trial_[i])) return invalid();
    staged_[i] = 1;
    return {};
  }
  Status stage_cphi(std::size_t cell, double value) noexcept {
    if (!preparing_ || cell >= cphi_.size() || inactive_[cell] ||
        !tcr::detail::dyn711_tick(clock()).update_cphi ||
        !std::isfinite(value) || value < 1 || value > 16) return invalid();
    trial_cphi_[cell] = value;
    cphi_staged_[cell] = 1;
    return {};
  }
  Status stage_inactive(std::size_t cell) noexcept {
    if (!preparing_ || cell >= cphi_.size()) return invalid();
    const auto begin = cell * ns_;
    std::copy(accepted_.begin()+begin, accepted_.begin()+begin+ns_, trial_.begin()+begin);
    trial_cphi_[cell] = cphi_[cell];
    std::fill(staged_.begin()+begin, staged_.begin()+begin+ns_, 1);
    cphi_staged_[cell] = inactive_[cell] = 1;
    return {};
  }
  Status seal() noexcept {
    if (!preparing_ || std::find(staged_.begin(), staged_.end(), 0) != staged_.end() ||
        std::find(cphi_staged_.begin(), cphi_staged_.end(), 0) != cphi_staged_.end() ||
        !valid(trial_, trial_cphi_)) return invalid();
    encode(trial_, trial_cphi_, candidate_step_, candidate_calls_, trial_bytes_);
    preparing_ = false;
    pending_ = true;
    return {};
  }
  void discard() noexcept { preparing_ = pending_ = false; }
  void commit() noexcept {
    if (!pending_) return;
    accepted_.swap(trial_); cphi_.swap(trial_cphi_); bytes_.swap(trial_bytes_);
    step_ = candidate_step_; calls_ = candidate_calls_;
    discard();
  }
  RestartCellRecordsView snapshot() const noexcept {
    return {identity_, width_, {bytes_.data(), bytes_.size()}};
  }
  RestartCellRecordsView prepared_snapshot() const noexcept {
    return pending_ ? RestartCellRecordsView{identity_, width_,
        {trial_bytes_.data(), trial_bytes_.size()}} : RestartCellRecordsView{};
  }
  Status restore(Span<const std::uint8_t> data, std::uint64_t step) noexcept {
    discard();
    if (data.size != bytes_.size() || !data.data || cphi_.empty()) return invalid();
    const auto calls = get(data.data + 8);
    if (calls > step) return invalid();
    const auto packed = pack(clock_at(calls));
    for (std::size_t cell=0; cell<cphi_.size(); ++cell) {
      const auto *p = data.data + cell * width_;
      if (get(p) != step || get(p+8) != calls || get(p+16) != ns_ || get(p+24) != packed)
        return invalid();
      trial_cphi_[cell] = real(p+32);
      for (std::size_t q=0; q<ns_; ++q) {
        const auto *r = p+40+40*q;
        const auto flags = get(r+32);
        if (flags > 2) return invalid(); // Upper and linear endpoint are exclusive.
        trial_[cell*ns_+q] = {real(r), real(r+8), real(r+16), real(r+24),
                             flags == 1, flags == 2};
      }
    }
    if (!valid(trial_, trial_cphi_)) return invalid();
    std::copy(data.data, data.data+data.size, trial_bytes_.begin());
    candidate_step_ = step; candidate_calls_ = calls; pending_ = true;
    return {};
  }
private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10241}; }
  static Clock clock_at(std::uint64_t calls) noexcept {
    return {static_cast<unsigned>(calls % 9),
        calls == 0 ? 0u : calls == 1 ? 6u : 12u-static_cast<unsigned>((calls-2)%7),
        calls >= 9};
  }
  static std::uint64_t pack(Clock c) noexcept {
    return c.rate_intervals | (std::uint64_t(c.cphi_count)<<8) |
        (std::uint64_t(c.rates_initialized)<<16);
  }
  static std::uint64_t get(const std::uint8_t *p) noexcept {
    std::uint64_t value{};
    for (unsigned i=0;i<8;++i) value |= std::uint64_t(p[i])<<(8*i);
    return value;
  }
  static void put(std::uint8_t *p, std::uint64_t value) noexcept {
    for (unsigned i=0;i<8;++i) p[i]=std::uint8_t(value>>(8*i));
  }
  static double real(const std::uint8_t *p) noexcept {
    const auto bits=get(p); double value; std::memcpy(&value,&bits,8); return value;
  }
  static void put_real(std::uint8_t *p, double value) noexcept {
    std::uint64_t bits; std::memcpy(&bits,&value,8); put(p,bits);
  }
  void encode(const std::vector<Rate> &values, const std::vector<double> &cphi,
      std::uint64_t step, std::uint64_t calls, std::vector<std::uint8_t> &bytes) const noexcept {
    for (std::size_t cell=0;cell<cphi.size();++cell) {
      auto *p=bytes.data()+cell*width_;
      put(p,step); put(p+8,calls); put(p+16,ns_); put(p+24,pack(clock_at(calls)));
      put_real(p+32,cphi[cell]);
      for (std::size_t q=0;q<ns_;++q) {
        auto *r=p+40+40*q; const auto &v=values[cell*ns_+q];
        put_real(r,v.pdf_sum); put_real(r+8,v.psr_sum);
        put_real(r+16,v.selected); put_real(r+24,v.effective);
        put(r+32,v.upper_branch ? 1 : v.linear_endpoint ? 2 : 0);
      }
    }
  }
  static bool valid(const std::vector<Rate> &values, const std::vector<double> &cphi) noexcept {
    for (const auto &v:values)
      if (!std::isfinite(v.pdf_sum) || !std::isfinite(v.psr_sum) ||
          !std::isfinite(v.selected) || v.selected < 0 ||
          !std::isfinite(v.effective) || v.effective != std::min(1.,v.selected) ||
          (v.upper_branch && v.linear_endpoint)) return false;
    for (double c:cphi) if (!std::isfinite(c) || c<1 || c>16) return false;
    return true;
  }
  PlanFingerprint identity_{};
  std::size_t ns_{};
  std::uint32_t width_{};
  std::uint64_t step_{}, calls_{}, candidate_step_{}, candidate_calls_{};
  bool preparing_{}, pending_{};
  std::vector<Rate> accepted_, trial_;
  std::vector<double> cphi_, trial_cphi_;
  std::vector<std::uint8_t> staged_, cphi_staged_, inactive_, bytes_, trial_bytes_;
};
} // namespace hundun::v04::detail

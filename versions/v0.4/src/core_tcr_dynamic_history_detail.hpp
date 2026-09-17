// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_io.hpp"
#include "models_tcr_dynamic_detail.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace hundun::v04::detail {
// Compact cell-major history. Slots per species: PDF/PSR interval rates,
// selected/effective kappa and discrete state. Three trailing slots hold Cd.
// Byte records follow global cells through the existing native Restart owner.
class DynamicTcrHistory {
public:
  void configure(PlanFingerprint model, std::size_t cells, std::size_t species,
                 double initial_cd, std::uint64_t initial_step = 0) {
    step_=initial_step; discard();
    ns_ = species;
    stride_ = 5 * ns_ + 3;
    width_ = static_cast<std::uint32_t>(24 + 8 * stride_);
    identity_ = (model ^ UINT64_C(0x54435244594e0001) ^ ns_) * UINT64_C(1099511628211);
    if (!identity_) identity_ = 1;
    accepted_.assign(cells * stride_, 0.);
    trial_.resize(accepted_.size());
    bytes_.resize(cells * width_);
    trial_bytes_.resize(bytes_.size());
    for (std::size_t i = 0; i < cells; ++i) {
      auto *p = accepted_.data() + i * stride_;
      for (std::size_t s = 0; s < ns_; ++s) p[5*s+2] = p[5*s+3] = 1.;
      std::fill(p + 5*ns_, p + stride_, initial_cd);
    }
    encode(accepted_, step_, bytes_);
  }
  std::uint64_t owned_bytes() const noexcept {
    return sizeof(*this) + 8 * (accepted_.capacity() + trial_.capacity()) +
           bytes_.capacity() + trial_bytes_.capacity();
  }
  std::size_t species() const noexcept { return ns_; }
  std::size_t stride() const noexcept { return stride_; }
  const double *accepted(std::size_t cell) const noexcept {
    return accepted_.data() + cell * stride_;
  }
  double *candidate(std::size_t cell) noexcept {
    return trial_.data() + cell * stride_;
  }
  Status begin(std::uint64_t step) noexcept {
    discard();
    if (step != step_ || step == UINT64_MAX) return invalid();
    std::copy(accepted_.begin(), accepted_.end(), trial_.begin());
    preparing_ = true;
    candidate_step_ = step + 1;
    return {};
  }
  Status seal() noexcept {
    if (!preparing_ || !valid(trial_)) return invalid();
    encode(trial_, candidate_step_, trial_bytes_);
    pending_ = true;
    preparing_ = false;
    return {};
  }
  void discard() noexcept { pending_ = preparing_ = false; }
  void commit() noexcept {
    if (!pending_) return;
    accepted_.swap(trial_);
    bytes_.swap(trial_bytes_);
    step_ = candidate_step_;
    discard();
  }
  RestartCellRecordsView snapshot() const noexcept {
    return {identity_, width_, {bytes_.data(), bytes_.size()}};
  }
  RestartCellRecordsView prepared_snapshot() const noexcept {
    return pending_ ? RestartCellRecordsView{identity_, width_,
                       {trial_bytes_.data(), trial_bytes_.size()}} : RestartCellRecordsView{};
  }
  Status restore(Span<const std::uint8_t> bytes, std::uint64_t step) noexcept {
    discard();
    if (bytes.size != bytes_.size() || (bytes.size && !bytes.data)) return invalid();
    for (std::size_t cell = 0; cell < accepted_.size()/stride_; ++cell) {
      const auto *p = bytes.data + cell * width_;
      if (get(p) != step || get(p+8) != step%4 || get(p+16) != ns_)
        return invalid();
      for (std::size_t i = 0; i < stride_; ++i) {
        const auto bits = get(p+24+8*i);
        std::memcpy(trial_.data()+cell*stride_+i, &bits, 8);
      }
    }
    if (!valid(trial_)) return invalid();
    std::copy(bytes.data, bytes.data+bytes.size, trial_bytes_.begin());
    candidate_step_ = step;
    pending_ = true;
    return {};
  }
private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10236}; }
  static std::uint64_t get(const std::uint8_t *p) noexcept {
    std::uint64_t v{}; for (unsigned i=0;i<8;++i) v|=std::uint64_t(p[i])<<(8*i);
    return v;
  }
  static void put(std::uint8_t *p, std::uint64_t v) noexcept {
    for (unsigned i=0;i<8;++i) p[i]=std::uint8_t(v>>(8*i));
  }
  void encode(const std::vector<double> &values, std::uint64_t step,
               std::vector<std::uint8_t> &bytes) const noexcept {
    for (std::size_t cell=0;cell<values.size()/stride_;++cell) {
      auto *p=bytes.data()+cell*width_;
      put(p,step); put(p+8,step%4); put(p+16,ns_);
      for (std::size_t i=0;i<stride_;++i) {
        std::uint64_t bits; std::memcpy(&bits,values.data()+cell*stride_+i,8);
        put(p+24+8*i,bits);
      }
    }
  }
  bool valid(const std::vector<double> &values) const noexcept {
    for (double v:values) if (!std::isfinite(v)) return false;
    for (std::size_t cell=0;cell<values.size()/stride_;++cell) {
      const auto *p=values.data()+cell*stride_;
      for (std::size_t s=0;s<ns_;++s)
        if (p[5*s+2]<0 || p[5*s+3]<1e-4 || p[5*s+3]>1 ||
            p[5*s+4]<0 || p[5*s+4]>3 || std::floor(p[5*s+4])!=p[5*s+4]) return false;
      for (unsigned g=0;g<3;++g)
        if (p[5*ns_+g]<1 || p[5*ns_+g]>16) return false;
    }
    return true;
  }
  PlanFingerprint identity_{};
  std::size_t ns_{}, stride_{};
  std::uint32_t width_{};
  std::uint64_t step_{}, candidate_step_{};
  bool preparing_{}, pending_{};
  std::vector<double> accepted_, trial_;
  std::vector<std::uint8_t> bytes_, trial_bytes_;
};
} // namespace hundun::v04::detail

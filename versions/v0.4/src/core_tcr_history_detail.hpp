// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_io.hpp"
#include "models_tcr_detail.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>

namespace hundun::v04::detail {
// Typed model history and its restart bytes share the native transaction.
// Staging can fail; publishing is only a pair of preallocated buffer swaps.
class ProductTcrHistory {
public:
  static constexpr std::uint32_t record_bytes = 120U;
  void configure(PlanFingerprint model, std::size_t cells,
                 int initialization_sign) {
    initialization_sign_ = initialization_sign;
    identity_ =
        (model ^ UINT64_C(0x5443524849530001)) * UINT64_C(1099511628211);
    if (identity_ == 0)
      identity_ = 1;
    accepted_.resize(cells);
    trial_.resize(cells);
    accepted_bytes_.resize(cells * record_bytes);
    trial_bytes_.resize(cells * record_bytes);
    for (std::size_t i = 0; i < cells; ++i) {
      accepted_[i].revision = {0, 1, 1};
      encode(accepted_[i], accepted_bytes_.data() + i * record_bytes);
    }
  }
  bool enabled() const noexcept { return identity_ != 0; }
  RestartCellRecordsView snapshot() const noexcept {
    return enabled() ? RestartCellRecordsView{identity_,
                                              record_bytes,
                                              {accepted_bytes_.data(),
                                               accepted_bytes_.size()}}
                     : RestartCellRecordsView{};
  }
  RestartCellRecordsView prepared_snapshot() const noexcept {
    return pending_ ? RestartCellRecordsView{identity_,
                                             record_bytes,
                                             {trial_bytes_.data(),
                                              trial_bytes_.size()}}
                    : RestartCellRecordsView{};
  }
  const tcr::detail::History &accepted(std::size_t cell) const noexcept {
    return accepted_[cell];
  }
  bool stage(std::size_t cell, const tcr::detail::Trial &trial,
             std::uint64_t step) noexcept {
    const auto &old = accepted_[cell];
    if (old.revision.accepted_step != step ||
        old.revision.input_revision == UINT64_MAX)
      return false;
    const auto candidate = tcr::detail::accept(
        old, trial, {step + 1, old.revision.input_revision + 1, 1});
    if (!candidate.available)
      return false;
    trial_[cell] = candidate.candidate;
    encode(trial_[cell], trial_bytes_.data() + cell * record_bytes);
    return true;
  }
  void seal() noexcept { pending_ = enabled(); }
  void discard() noexcept { pending_ = false; }
  void commit() noexcept {
    assert(!enabled() || pending_);
    if (!pending_)
      return;
    accepted_.swap(trial_);
    accepted_bytes_.swap(trial_bytes_);
    pending_ = false;
  }
  Status stage_restore(const RestartImage &image) noexcept {
    discard();
    if (!image.cell_record_lengths.empty() ||
        image.cell_record_identity != identity_ ||
        image.cell_record_bytes != (enabled() ? record_bytes : 0U) ||
        image.cell_records.size() != accepted_bytes_.size())
      return invalid();
    if ((image.source_format_version == 4) != enabled())
      return invalid();
    if (!enabled())
      return {};
    if (image.source_format_version != 4 || image.backward_euler_recovery)
      return invalid();
    return stage_restore_records(
        {image.cell_records.data(), image.cell_records.size()}, image.step);
  }
  // A combined V5 owner has already checked its outer model identity/format
  // and exact native clock; this still validates every typed TCR history.
  Status stage_restore_records(Span<const std::uint8_t> records,
                               std::uint64_t step) noexcept {
    discard();
    if (records.size != accepted_bytes_.size() ||
        (records.size && !records.data))
      return invalid();
    if (!enabled())
      return {};
    for (std::size_t i = 0; i < accepted_.size(); ++i) {
      auto &history = trial_[i];
      if (!decode(records.data + i * record_bytes, history) ||
          history.revision.accepted_step != step ||
          history.revision.input_revision == 0 ||
          !tcr::detail::restore(history, history.revision).available ||
          (history.initialized &&
           (history.mapping_identity !=
                tcr::detail::kReactantMoleFractionMappingIdentity ||
            (initialization_sign_ != 0 &&
             history.initialization_sign != initialization_sign_))))
        return invalid();
    }
    std::copy(records.data, records.data + records.size, trial_bytes_.begin());
    seal();
    return {};
  }

private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10217}; }
  static void put(std::uint8_t *&p, std::uint64_t v, unsigned n) noexcept {
    for (unsigned i = 0; i < n; ++i)
      *p++ = std::uint8_t(v >> (8 * i));
  }
  static std::uint64_t get(const std::uint8_t *&p, unsigned n) noexcept {
    std::uint64_t value = 0;
    for (unsigned i = 0; i < n; ++i)
      value |= std::uint64_t(*p++) << (8 * i);
    return value;
  }
  static void real(std::uint8_t *&p, double v) noexcept {
    std::uint64_t bits;
    std::memcpy(&bits, &v, 8);
    put(p, bits, 8);
  }
  static double real(const std::uint8_t *&p) noexcept {
    const auto bits = get(p, 8);
    double value;
    std::memcpy(&value, &bits, 8);
    return value;
  }
  static void revision(std::uint8_t *&p, portable::Revision v) noexcept {
    put(p, v.accepted_step, 8);
    put(p, v.input_revision, 8);
    put(p, v.algorithm_version, 4);
  }
  static portable::Revision revision(const std::uint8_t *&p) noexcept {
    const auto step = get(p, 8), input = get(p, 8);
    const auto algorithm = get(p, 4);
    return {step, input, std::uint32_t(algorithm)};
  }
  static void encode(const tcr::detail::History &h, std::uint8_t *p) noexcept {
    revision(p, h.revision);
    put(p, h.initialized, 4);
    real(p, h.input.eta);
    real(p, h.input.rate_ratio);
    real(p, h.signed_root);
    real(p, h.control);
    put(p, h.branch_sign + 1, 4);
    put(p, h.initialization_sign + 1, 4);
    put(p, h.mapping_identity, 8);
    put(p, h.fold_count, 8);
    real(p, h.last_fold.eta);
    real(p, h.last_fold.rate_ratio);
    revision(p, h.last_fold_base_revision);
    put(p, 0, 4);
  }
  static bool decode(const std::uint8_t *p, tcr::detail::History &h) noexcept {
    h = {};
    h.revision = revision(p);
    const auto initialized = get(p, 4);
    h.initialized = initialized != 0;
    h.input = {real(p), real(p)};
    h.signed_root = real(p);
    h.control = real(p);
    const auto branch = get(p, 4), initial = get(p, 4);
    if (initialized > 1 || branch > 2 || initial > 2)
      return false;
    h.branch_sign = int(branch) - 1;
    h.initialization_sign = int(initial) - 1;
    h.mapping_identity = get(p, 8);
    h.fold_count = get(p, 8);
    h.last_fold = {real(p), real(p)};
    h.last_fold_base_revision = revision(p);
    if (get(p, 4) != 0)
      return false;
    // A nonexistent fold has one canonical encoding, including its revision.
    if (h.fold_count == 0 &&
        (h.last_fold.eta != 0 || h.last_fold.rate_ratio != 0 ||
         h.last_fold_base_revision != portable::Revision{}))
      return false;
    return true;
  }
  int initialization_sign_{};
  PlanFingerprint identity_{};
  bool pending_{};
  std::vector<tcr::detail::History> accepted_, trial_;
  std::vector<std::uint8_t> accepted_bytes_, trial_bytes_;
};
} // namespace hundun::v04::detail

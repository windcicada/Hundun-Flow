// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_io.hpp"
#include "models_spray_migration_detail.hpp"
#include "models_spray_parcel_detail.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

namespace hundun::v04::detail {
// Native variable-cell records own parcel/TAB/lineage values and borrow the
// product's configured injectors. A pending record set, parcel set and injector
// counter set are published only by the native gas transaction's final commit.
class ProductSprayHistory {
public:
  using Parcel = spray::detail::ParcelMigrationValue;
  using Injector = spray::detail::DeterministicInjector;
  static constexpr std::size_t header_bytes = 24, parcel_bytes = 144,
                               injector_bytes = 24;
  Status configure(PlanFingerprint identity, MeshPatch patch, Int3 global,
                   std::size_t capacity, std::uint64_t material,
                   Span<Injector *const> injectors, RestartCellRecordsView tcr,
                   std::uint64_t maximum_bytes) {
    if (identity_ || !identity || !material || !capacity || global.x <= 0 ||
        global.y <= 0 || global.z <= 0 || patch.cells.x <= 0 ||
        patch.cells.y <= 0 || patch.cells.z <= 0 || injectors.size > 64 ||
        (injectors.size && !injectors.data) ||
        (tcr.identity ? tcr.record_bytes != 120 : tcr.record_bytes != 0))
      return invalid();
    global_count_ = 1;
    count_ = 1;
    const int extents[]{global.x, global.y, global.z};
    const int local[]{patch.cells.x, patch.cells.y, patch.cells.z};
    const int begin[]{patch.begin.x, patch.begin.y, patch.begin.z};
    for (unsigned d = 0; d < 3; ++d) {
      if (begin[d] < 0 || local[d] > extents[d] ||
          begin[d] > extents[d] - local[d] ||
          global_count_ > UINT64_MAX / std::uint64_t(extents[d]) ||
          count_ > SIZE_MAX / std::size_t(local[d]))
        return invalid();
      global_count_ *= extents[d];
      count_ *= local[d];
    }
    global_ = global;
    patch_ = patch;
    capacity_ = capacity;
    material_ = material;
    tcr_identity_ = tcr.identity;
    tcr_width_ = tcr.record_bytes;
    if (capacity > (UINT32_MAX - header_bytes - tcr_width_) / parcel_bytes ||
        count_ > (SIZE_MAX - capacity * parcel_bytes -
                  injectors.size * injector_bytes) /
                     (header_bytes + tcr_width_))
      return invalid();
    maximum_payload_ = count_ * (header_bytes + tcr_width_) +
                       capacity * parcel_bytes +
                       injectors.size * injector_bytes;
    // Conservative resident bound for both typed/encoded state sets and
    // indexing scratch. Borrowed injector storage is accounted by its owner.
    const auto fixed =
        sizeof(*this) +
        injectors.size * (sizeof(Injector *) +
                          sizeof(spray::detail::InjectorCommittedState) +
                          sizeof(std::uint8_t));
    if (fixed > maximum_bytes)
      return {StatusCode::allocation_failure, 10241};
    const auto supplied_budget = maximum_bytes;
    maximum_bytes -= fixed;
    if (maximum_payload_ > maximum_bytes / 2 ||
        capacity >
            (maximum_bytes - 2 * maximum_payload_) / (2 * sizeof(Parcel)) ||
        count_ >
            (maximum_bytes - 2 * maximum_payload_ -
             2 * capacity * sizeof(Parcel)) /
                (2 * sizeof(std::uint32_t) + sizeof(std::size_t) + tcr_width_))
      return {StatusCode::allocation_failure, 10241};
    owned_bytes_ =
        fixed + 2 * maximum_payload_ + 2 * capacity * sizeof(Parcel) +
        count_ * (2 * sizeof(std::uint32_t) + sizeof(std::size_t) + tcr_width_);
    if (owned_bytes_ > supplied_budget)
      return invalid();
    if (injectors.size)
      injectors_.assign(injectors.data, injectors.data + injectors.size);
    for (std::size_t i = 0; i < injectors_.size(); ++i) {
      std::size_t local{};
      if (!injectors_[i] || !injectors_[i]->configured() ||
          injectors_[i]->trial_active() ||
          injectors_[i]->configured_spec().liquid_material_fingerprint !=
              material_ ||
          !local_cell(injectors_[i]->configured_spec().owner_global_cell,
                      local))
        return invalid();
      for (std::size_t j = 0; j < i; ++j)
        if (injectors_[j]->configured_spec().injector_id ==
            injectors_[i]->configured_spec().injector_id)
          return invalid();
    }
    accepted_.reserve(capacity);
    trial_.reserve(capacity);
    accepted_bytes_.reserve(maximum_payload_);
    trial_bytes_.reserve(maximum_payload_);
    accepted_lengths_.resize(count_);
    trial_lengths_.resize(count_);
    offsets_.resize(count_);
    trial_tcr_.resize(count_ * tcr_width_);
    injector_states_.resize(injectors_.size());
    injector_seen_.resize(injectors_.size());
    if (!encode(accepted_, 0, tcr, false, accepted_lengths_, accepted_bytes_))
      return invalid();
    identity_ = identity;
    return {};
  }
  std::uint64_t owned_bytes() const noexcept { return owned_bytes_; }
  std::size_t maximum_payload_bytes() const noexcept {
    return maximum_payload_;
  }
  bool enabled() const noexcept { return identity_ != 0; }
  Span<const Parcel> accepted_parcels() const noexcept {
    return {accepted_.data(), accepted_.size()};
  }
  Span<const Parcel> prepared_parcels() const noexcept {
    return pending_ ? Span<const Parcel>{trial_.data(), trial_.size()}
                    : Span<const Parcel>{};
  }
  Span<const std::uint8_t> prepared_tcr_records() const noexcept {
    return pending_
               ? Span<const std::uint8_t>{trial_tcr_.data(), trial_tcr_.size()}
               : Span<const std::uint8_t>{};
  }
  RestartCellRecordsView snapshot() const noexcept {
    return {identity_,
            0,
            {accepted_bytes_.data(), accepted_bytes_.size()},
            {accepted_lengths_.data(), accepted_lengths_.size()}};
  }
  Status stage_next(Span<const Parcel> parcels, std::uint64_t accepted_step,
                    RestartCellRecordsView tcr) noexcept {
    if (!enabled() || pending_ || accepted_step != step_ ||
        step_ == UINT64_MAX || parcels.size > capacity_ ||
        (parcels.size && !parcels.data))
      return invalid();
    trial_.clear();
    if (parcels.size)
      trial_.assign(parcels.data, parcels.data + parcels.size);
    if (!valid_parcels() ||
        !encode(trial_, step_ + 1, tcr, true, trial_lengths_, trial_bytes_)) {
      discard();
      return invalid();
    }
    if (tcr.values.size)
      std::copy(tcr.values.data, tcr.values.data + tcr.values.size,
                trial_tcr_.begin());
    pending_step_ = step_ + 1;
    pending_ = true;
    return {};
  }
  Status stage_restore(const RestartImage &image) noexcept {
    discard();
    if (!enabled() || image.source_format_version != 5 ||
        image.backward_euler_recovery ||
        image.cell_record_identity != identity_ ||
        image.cell_record_bytes != 0 ||
        image.cell_record_lengths.size() != count_ ||
        image.cell_records.size() > maximum_payload_)
      return invalid();
    trial_.clear();
    std::fill(injector_seen_.begin(), injector_seen_.end(), false);
    std::size_t offset = 0;
    for (std::size_t cell = 0; cell < count_; ++cell) {
      const std::size_t length = image.cell_record_lengths[cell];
      if (offset > image.cell_records.size() ||
          length > image.cell_records.size() - offset)
        return invalid();
      if (length == 0) {
        if (tcr_width_)
          return invalid();
        continue;
      }
      if (length < header_bytes)
        return invalid();
      const auto *p = image.cell_records.data() + offset;
      if (get(p) != image.step)
        return invalid();
      const auto np = get(p, 4), ni = get(p, 4), nt = get(p, 4),
                 version = get(p, 4);
      if (version != 1 || nt != tcr_width_ || np > capacity_ - trial_.size() ||
          ni > injectors_.size() ||
          header_bytes + nt + np * parcel_bytes + ni * injector_bytes != length)
        return invalid();
      if (nt) {
        std::memcpy(trial_tcr_.data() + cell * tcr_width_, p, nt);
        p += nt;
      }
      for (std::uint64_t i = 0; i < np; ++i) {
        const auto value = decode_parcel(p);
        std::size_t local{};
        if (!local_cell(value.parcel.owner_global_cell, local) || local != cell)
          return invalid();
        trial_.push_back(value);
      }
      for (std::uint64_t i = 0; i < ni; ++i) {
        const auto id = get(p);
        spray::detail::InjectorCommittedState value{real(p), get(p)};
        std::size_t index = 0;
        for (; index < injectors_.size(); ++index)
          if (injectors_[index]->configured_spec().injector_id == id)
            break;
        std::size_t local{};
        if (index == injectors_.size() || injector_seen_[index] ||
            !local_cell(injectors_[index]->configured_spec().owner_global_cell,
                        local) ||
            local != cell)
          return invalid();
        injector_states_[index] = value;
        injector_seen_[index] = true;
      }
      offset += length;
    }
    if (offset != image.cell_records.size() ||
        std::find(injector_seen_.begin(), injector_seen_.end(), false) !=
            injector_seen_.end() ||
        !valid_parcels())
      return invalid();
    for (std::size_t i = 0; i < injectors_.size(); ++i)
      if (!injectors_[i]->stage_restore(injector_states_[i])) {
        discard();
        return invalid();
      }
    trial_bytes_.assign(image.cell_records.begin(), image.cell_records.end());
    std::copy(image.cell_record_lengths.begin(),
              image.cell_record_lengths.end(), trial_lengths_.begin());
    pending_step_ = image.step;
    pending_ = true;
    return {};
  }
  Status preflight_commit() const noexcept {
    if (!pending_)
      return invalid();
    for (auto *injector : injectors_)
      if (!injector->preflight_commit())
        return invalid();
    return {};
  }
  void commit() noexcept {
    assert(preflight_commit());
    for (auto *injector : injectors_)
      injector->publish_preflighted_trial();
    accepted_.swap(trial_);
    accepted_bytes_.swap(trial_bytes_);
    accepted_lengths_.swap(trial_lengths_);
    step_ = pending_step_;
    pending_ = false;
  }
  void discard() noexcept {
    for (auto *injector : injectors_)
      if (injector && injector->trial_active())
        (void)injector->rollback_trial();
    pending_ = false;
  }

private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10241}; }
  bool local_cell(std::uint64_t id, std::size_t &out) const noexcept {
    const std::uint64_t nx = global_.x, ny = global_.y;
    if (id >= global_count_)
      return false;
    const auto x = std::int64_t(id % nx) - patch_.begin.x;
    const auto y = std::int64_t((id / nx) % ny) - patch_.begin.y;
    const auto z = std::int64_t(id / (nx * ny)) - patch_.begin.z;
    if (x < 0 || y < 0 || z < 0 || x >= patch_.cells.x || y >= patch_.cells.y ||
        z >= patch_.cells.z)
      return false;
    out = std::size_t(x) +
          patch_.cells.x * (std::size_t(y) + patch_.cells.y * std::size_t(z));
    return true;
  }
  bool valid_parcels() noexcept {
    for (const auto &v : trial_) {
      std::size_t cell{};
      if (spray::validate_parcel_state(v.parcel) !=
              spray::ParcelStateStatus::success ||
          v.parcel.liquid_material_fingerprint != material_ ||
          !std::isfinite(v.tab_deformation) ||
          !std::isfinite(v.tab_deformation_rate_per_s) ||
          !local_cell(v.parcel.owner_global_cell, cell))
        return false;
    }
    std::sort(trial_.begin(), trial_.end(),
              [](const Parcel &a, const Parcel &b) {
                return a.parcel.id < b.parcel.id;
              });
    for (std::size_t i = 1; i < trial_.size(); ++i)
      if (trial_[i - 1].parcel.id == trial_[i].parcel.id)
        return false;
    std::sort(
        trial_.begin(), trial_.end(), [](const Parcel &a, const Parcel &b) {
          return a.parcel.owner_global_cell < b.parcel.owner_global_cell ||
                 (a.parcel.owner_global_cell == b.parcel.owner_global_cell &&
                  a.parcel.id < b.parcel.id);
        });
    return true;
  }
  bool encode(const std::vector<Parcel> &parcels, std::uint64_t step,
              RestartCellRecordsView tcr, bool pending_injectors,
              std::vector<std::uint32_t> &lengths,
              std::vector<std::uint8_t> &bytes) noexcept {
    if (tcr.identity != tcr_identity_ || tcr.record_bytes != tcr_width_ ||
        tcr.values.size != count_ * tcr_width_ ||
        (tcr.values.size && !tcr.values.data) || tcr.variable_cell_bytes.size)
      return false;
    std::fill(lengths.begin(), lengths.end(),
              tcr_width_ ? header_bytes + tcr_width_ : 0);
    const auto add = [&](std::uint64_t global, std::size_t amount) {
      std::size_t cell{};
      if (!local_cell(global, cell))
        return false;
      auto &n = lengths[cell];
      if (!n)
        n = header_bytes;
      if (amount > UINT32_MAX - n)
        return false;
      n += std::uint32_t(amount);
      return true;
    };
    for (const auto &p : parcels)
      if (!add(p.parcel.owner_global_cell, parcel_bytes))
        return false;
    for (auto *i : injectors_)
      if (!add(i->configured_spec().owner_global_cell, injector_bytes))
        return false;
    std::size_t total = 0;
    for (std::size_t cell = 0; cell < count_; ++cell) {
      offsets_[cell] = total;
      if (lengths[cell] > maximum_payload_ - total)
        return false;
      total += lengths[cell];
    }
    bytes.resize(total);
    for (std::size_t cell = 0; cell < count_; ++cell) {
      if (!lengths[cell])
        continue;
      auto *p = bytes.data() + offsets_[cell];
      put(p, step);
      put(p, 0, 4);
      put(p, 0, 4);
      put(p, tcr_width_, 4);
      put(p, 1, 4);
      if (tcr_width_) {
        const auto *source = tcr.values.data + cell * tcr_width_;
        const auto *cursor = source;
        if (get(cursor) != step)
          return false;
        std::memcpy(p, source, tcr_width_);
        p += tcr_width_;
      }
      offsets_[cell] = std::size_t(p - bytes.data());
    }
    for (const auto &v : parcels) {
      std::size_t cell{};
      if (!local_cell(v.parcel.owner_global_cell, cell))
        return false;
      auto *p = bytes.data() + offsets_[cell];
      encode_parcel(v, p);
      offsets_[cell] += parcel_bytes;
    }
    for (auto *i : injectors_) {
      spray::detail::InjectorCommittedState state = i->committed_state();
      if (pending_injectors && !i->prepared_state(state))
        return false;
      std::size_t cell{};
      if (!local_cell(i->configured_spec().owner_global_cell, cell))
        return false;
      auto *p = bytes.data() + offsets_[cell];
      put(p, i->configured_spec().injector_id);
      real(p, state.residual_mass_kg);
      put(p, state.next_ordinal);
      offsets_[cell] += injector_bytes;
    }
    // Counts are derived from the final byte extents, with one linear pass
    // over the cell-sorted parcels and the bounded injector set.
    std::size_t offset = 0, parcel = 0;
    for (std::size_t cell = 0; cell < count_; ++cell) {
      if (!lengths[cell])
        continue;
      std::uint32_t np = 0, ni = 0;
      while (parcel < parcels.size()) {
        std::size_t owner{};
        if (!local_cell(parcels[parcel].parcel.owner_global_cell, owner))
          return false;
        if (owner != cell)
          break;
        ++np;
        ++parcel;
      }
      for (auto *i : injectors_) {
        std::size_t owner{};
        if (!local_cell(i->configured_spec().owner_global_cell, owner))
          return false;
        ni += owner == cell;
      }
      auto *p = bytes.data() + offset + 8;
      put(p, np, 4);
      put(p, ni, 4);
      offset += lengths[cell];
    }
    return parcel == parcels.size();
  }
  static void put(std::uint8_t *&p, std::uint64_t value,
                  unsigned n = 8) noexcept {
    for (unsigned i = 0; i < n; ++i)
      *p++ = std::uint8_t(value >> (8 * i));
  }
  static std::uint64_t get(const std::uint8_t *&p, unsigned n = 8) noexcept {
    std::uint64_t value = 0;
    for (unsigned i = 0; i < n; ++i)
      value |= std::uint64_t(*p++) << (8 * i);
    return value;
  }
  static void real(std::uint8_t *&p, double value) noexcept {
    std::uint64_t bits;
    std::memcpy(&bits, &value, 8);
    put(p, bits);
  }
  static double real(const std::uint8_t *&p) noexcept {
    const auto bits = get(p);
    double value;
    std::memcpy(&value, &bits, 8);
    return value;
  }
  static void encode_parcel(const Parcel &v, std::uint8_t *&p) noexcept {
    const auto &q = v.parcel;
    put(p, q.id.high);
    put(p, q.id.low);
    for (double x : q.position_m)
      real(p, x);
    for (double x : q.velocity_m_per_s)
      real(p, x);
    real(p, q.droplet_mass_kg);
    real(p, q.droplet_diameter_m);
    real(p, q.multiplicity);
    real(p, q.temperature_k);
    put(p, q.liquid_material_fingerprint);
    put(p, q.owner_global_cell);
    real(p, q.age_s);
    real(p, v.tab_deformation);
    real(p, v.tab_deformation_rate_per_s);
    put(p, v.breakup_ordinal);
  }
  static Parcel decode_parcel(const std::uint8_t *&p) noexcept {
    Parcel v;
    auto &q = v.parcel;
    q.id = {get(p), get(p)};
    for (double &x : q.position_m)
      x = real(p);
    for (double &x : q.velocity_m_per_s)
      x = real(p);
    q.droplet_mass_kg = real(p);
    q.droplet_diameter_m = real(p);
    q.multiplicity = real(p);
    q.temperature_k = real(p);
    q.liquid_material_fingerprint = get(p);
    q.owner_global_cell = get(p);
    q.age_s = real(p);
    v.tab_deformation = real(p);
    v.tab_deformation_rate_per_s = real(p);
    v.breakup_ordinal = get(p);
    return v;
  }
  PlanFingerprint identity_{}, tcr_identity_{};
  MeshPatch patch_{};
  Int3 global_{};
  std::size_t count_{}, capacity_{}, maximum_payload_{};
  std::uint64_t material_{}, step_{}, pending_step_{}, global_count_{},
      owned_bytes_{};
  std::uint32_t tcr_width_{};
  std::vector<Injector *> injectors_;
  std::vector<spray::detail::InjectorCommittedState> injector_states_;
  std::vector<std::uint8_t> injector_seen_;
  std::vector<Parcel> accepted_, trial_;
  std::vector<std::uint8_t> accepted_bytes_, trial_bytes_, trial_tcr_;
  std::vector<std::uint32_t> accepted_lengths_, trial_lengths_;
  std::vector<std::size_t> offsets_;
  bool pending_{};
};
} // namespace hundun::v04::detail

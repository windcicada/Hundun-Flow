// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_spray.hpp"
#include "hundun/v04_status.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace hundun::v04::spray::detail {

enum class ParcelOperationDetail : std::uint32_t {
  none = 0U,
  trial_already_active = 1U,
  trial_not_active = 2U,
  invalid_parcel = 3U,
  duplicate_id = 4U,
  capacity_exceeded = 5U,
  parcel_not_found = 6U,
  invalid_snapshot = 7U,
  invalid_injector = 8U,
  injector_not_configured = 9U,
  injector_counter_overflow = 10U,
  invalid_tab_state = 11U
};

struct ParcelSoASnapshot {
  std::vector<std::uint64_t> id_high;
  std::vector<std::uint64_t> id_low;
  std::vector<double> position_x_m;
  std::vector<double> position_y_m;
  std::vector<double> position_z_m;
  std::vector<double> velocity_x_m_per_s;
  std::vector<double> velocity_y_m_per_s;
  std::vector<double> velocity_z_m_per_s;
  std::vector<double> droplet_mass_kg;
  std::vector<double> droplet_diameter_m;
  std::vector<double> multiplicity;
  std::vector<double> temperature_k;
  std::vector<std::uint64_t> liquid_material_fingerprint;
  std::vector<std::uint64_t> owner_global_cell;
  std::vector<double> age_s;
  std::vector<double> tab_deformation;
  std::vector<double> tab_deformation_rate_per_s;

  [[nodiscard]] std::size_t size() const noexcept { return id_high.size(); }
  [[nodiscard]] bool consistent() const noexcept;
};

class ParcelContainer {
public:
  ParcelContainer() = default;
  ParcelContainer(const ParcelContainer&) = delete;
  ParcelContainer& operator=(const ParcelContainer&) = delete;
  ParcelContainer(ParcelContainer&&) = delete;
  ParcelContainer& operator=(ParcelContainer&&) = delete;

  [[nodiscard]] Status reserve(std::size_t maximum_parcels) noexcept;

  [[nodiscard]] std::size_t capacity() const noexcept {
    return capacity_limit_;
  }
  [[nodiscard]] std::size_t committed_size() const noexcept {
    return committed_.size();
  }
  [[nodiscard]] std::size_t trial_size() const noexcept {
    return trial_active_ ? trial_.size() : 0U;
  }
  [[nodiscard]] bool trial_active() const noexcept { return trial_active_; }

  [[nodiscard]] bool committed_at(std::size_t index,
                                  SprayParcelState& out) const noexcept;
  [[nodiscard]] bool trial_at(std::size_t index,
                              SprayParcelState& out) const noexcept;
  [[nodiscard]] bool committed_tab_state(
      ParcelId id, double& deformation,
      double& deformation_rate_per_s) const noexcept;
  [[nodiscard]] bool trial_tab_state(
      ParcelId id, double& deformation,
      double& deformation_rate_per_s) const noexcept;

  [[nodiscard]] Status begin_trial() noexcept;
  [[nodiscard]] Status stage_add(const SprayParcelState& parcel) noexcept;
  [[nodiscard]] Status stage_update(const SprayParcelState& parcel) noexcept;
  [[nodiscard]] Status stage_tab_state(
      ParcelId id, double deformation,
      double deformation_rate_per_s) noexcept;
  [[nodiscard]] Status stage_remove(ParcelId id) noexcept;
  [[nodiscard]] Status preflight_commit() const noexcept;
  // Precondition: preflight_commit() just succeeded and no intervening
  // mutation occurred. A synchronized shared-finish wrapper normally calls
  // this no-fail publication step after every participant preflights.
  void publish_preflighted_trial() noexcept;
  [[nodiscard]] Status commit_trial() noexcept;
  [[nodiscard]] Status rollback_trial() noexcept;

  [[nodiscard]] Status
  snapshot_committed(ParcelSoASnapshot& out) const noexcept;
  [[nodiscard]] Status
  restore_committed(const ParcelSoASnapshot& snapshot) noexcept;

private:
  [[nodiscard]] Status latch(Status failure) noexcept;

  ParcelSoASnapshot committed_;
  ParcelSoASnapshot trial_;
  std::size_t capacity_limit_{};
  bool trial_active_{};
  Status trial_failure_{};
};

enum class InjectorShape : std::uint8_t { point, cone };

struct InjectorSpec {
  std::uint64_t seed{};
  std::uint64_t injector_id{};
  InjectorShape shape{InjectorShape::point};
  Vector3 origin_m{};
  Vector3 axis{1.0, 0.0, 0.0};
  double cone_half_angle_rad{};
  double injection_speed_m_per_s{};
  double mass_flow_rate_kg_per_s{};
  double represented_mass_per_parcel_kg{};
  double droplet_mass_kg{};
  double droplet_diameter_m{};
  double temperature_k{};
  std::uint64_t liquid_material_fingerprint{};
  std::uint64_t owner_global_cell{};
};

struct InjectorCommittedState {
  double residual_mass_kg{};
  std::uint64_t next_ordinal{};

  friend constexpr bool operator==(InjectorCommittedState left,
                                   InjectorCommittedState right) noexcept {
    return left.residual_mass_kg == right.residual_mass_kg &&
           left.next_ordinal == right.next_ordinal;
  }
  friend constexpr bool operator!=(InjectorCommittedState left,
                                   InjectorCommittedState right) noexcept {
    return !(left == right);
  }
};

enum class InjectionStatus : std::uint8_t {
  success,
  invalid_input,
  not_configured,
  trial_already_active,
  capacity_exceeded,
  counter_overflow,
  non_finite_output
};

struct InjectionReport {
  InjectionStatus status{InjectionStatus::invalid_input};
  std::string_view model_id{"deterministic_point_cone_injector_v1"};
  std::uint64_t accepted_step{};
  std::uint64_t injector_id{};
  std::uint64_t first_ordinal{};
  std::uint64_t next_ordinal{};
  std::size_t parcel_count{};
  double requested_mass_kg{};
  double injected_mass_kg{};
  double residual_mass_before_kg{};
  double residual_mass_after_kg{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == InjectionStatus::success;
  }
};

class DeterministicInjector {
public:
  DeterministicInjector() = default;
  DeterministicInjector(const DeterministicInjector&) = delete;
  DeterministicInjector& operator=(const DeterministicInjector&) = delete;
  DeterministicInjector(DeterministicInjector&&) = delete;
  DeterministicInjector& operator=(DeterministicInjector&&) = delete;

  [[nodiscard]] Status reserve(std::size_t maximum_candidates) noexcept;
  [[nodiscard]] Status configure(
      const InjectorSpec& spec,
      InjectorCommittedState initial_state = {}) noexcept;

  [[nodiscard]] InjectionReport begin_trial(std::uint64_t accepted_step,
                                            double duration_s) noexcept;
  [[nodiscard]] Status preflight_commit() const noexcept;
  // Precondition: preflight_commit() just succeeded and no intervening
  // mutation occurred. A synchronized shared-finish wrapper normally calls
  // this no-fail publication step after every participant preflights.
  void publish_preflighted_trial() noexcept;
  [[nodiscard]] Status commit_trial() noexcept;
  [[nodiscard]] Status rollback_trial() noexcept;

  [[nodiscard]] bool configured() const noexcept { return configured_; }
  [[nodiscard]] bool trial_active() const noexcept { return trial_active_; }
  [[nodiscard]] std::size_t candidate_capacity() const noexcept {
    return candidate_capacity_;
  }
  [[nodiscard]] std::size_t candidate_count() const noexcept {
    return trial_active_ ? candidates_.size() : 0U;
  }
  [[nodiscard]] bool candidate_at(std::size_t index,
                                  SprayParcelState& out) const noexcept;
  [[nodiscard]] InjectorCommittedState committed_state() const noexcept {
    return committed_;
  }
  [[nodiscard]] const InjectionReport& trial_report() const noexcept {
    return report_;
  }

private:
  InjectorSpec spec_{};
  InjectorCommittedState committed_{};
  InjectorCommittedState trial_state_{};
  ParcelSoASnapshot candidates_;
  std::vector<ParcelId> candidate_id_scratch_;
  std::size_t candidate_capacity_{};
  bool configured_{};
  bool trial_active_{};
  InjectionReport report_{};
};

} // namespace hundun::v04::spray::detail

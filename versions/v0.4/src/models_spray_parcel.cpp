// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_parcel_detail.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMaximumExactlyRepresentableCount = 9007199254740991.0;

Status failure(StatusCode code, ParcelOperationDetail detail) noexcept {
  return {code, static_cast<std::uint32_t>(detail)};
}

bool finite_vector(const Vector3& value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}

bool valid_id(ParcelId id) noexcept {
  return id.high != 0U || id.low != 0U;
}

ParcelId id_at(const ParcelSoASnapshot& arrays, std::size_t index) noexcept {
  return {arrays.id_high[index], arrays.id_low[index]};
}

bool arrays_have_capacity(const ParcelSoASnapshot& arrays,
                          std::size_t count) noexcept {
  return arrays.id_high.capacity() >= count &&
         arrays.id_low.capacity() >= count &&
         arrays.position_x_m.capacity() >= count &&
         arrays.position_y_m.capacity() >= count &&
         arrays.position_z_m.capacity() >= count &&
         arrays.velocity_x_m_per_s.capacity() >= count &&
         arrays.velocity_y_m_per_s.capacity() >= count &&
         arrays.velocity_z_m_per_s.capacity() >= count &&
         arrays.droplet_mass_kg.capacity() >= count &&
         arrays.droplet_diameter_m.capacity() >= count &&
         arrays.multiplicity.capacity() >= count &&
         arrays.temperature_k.capacity() >= count &&
         arrays.liquid_material_fingerprint.capacity() >= count &&
         arrays.owner_global_cell.capacity() >= count &&
         arrays.age_s.capacity() >= count &&
         arrays.tab_deformation.capacity() >= count &&
         arrays.tab_deformation_rate_per_s.capacity() >= count;
}

void reserve_arrays(ParcelSoASnapshot& arrays, std::size_t count) {
  arrays.id_high.reserve(count);
  arrays.id_low.reserve(count);
  arrays.position_x_m.reserve(count);
  arrays.position_y_m.reserve(count);
  arrays.position_z_m.reserve(count);
  arrays.velocity_x_m_per_s.reserve(count);
  arrays.velocity_y_m_per_s.reserve(count);
  arrays.velocity_z_m_per_s.reserve(count);
  arrays.droplet_mass_kg.reserve(count);
  arrays.droplet_diameter_m.reserve(count);
  arrays.multiplicity.reserve(count);
  arrays.temperature_k.reserve(count);
  arrays.liquid_material_fingerprint.reserve(count);
  arrays.owner_global_cell.reserve(count);
  arrays.age_s.reserve(count);
  arrays.tab_deformation.reserve(count);
  arrays.tab_deformation_rate_per_s.reserve(count);
}

void resize_arrays(ParcelSoASnapshot& arrays, std::size_t count) noexcept {
  arrays.id_high.resize(count);
  arrays.id_low.resize(count);
  arrays.position_x_m.resize(count);
  arrays.position_y_m.resize(count);
  arrays.position_z_m.resize(count);
  arrays.velocity_x_m_per_s.resize(count);
  arrays.velocity_y_m_per_s.resize(count);
  arrays.velocity_z_m_per_s.resize(count);
  arrays.droplet_mass_kg.resize(count);
  arrays.droplet_diameter_m.resize(count);
  arrays.multiplicity.resize(count);
  arrays.temperature_k.resize(count);
  arrays.liquid_material_fingerprint.resize(count);
  arrays.owner_global_cell.resize(count);
  arrays.age_s.resize(count);
  arrays.tab_deformation.resize(count);
  arrays.tab_deformation_rate_per_s.resize(count);
}

template <class T>
void copy_vector_noalloc(std::vector<T>& destination,
                         const std::vector<T>& source) noexcept {
  destination.resize(source.size());
  std::copy(source.begin(), source.end(), destination.begin());
}

void copy_arrays_noalloc(ParcelSoASnapshot& destination,
                         const ParcelSoASnapshot& source) noexcept {
  copy_vector_noalloc(destination.id_high, source.id_high);
  copy_vector_noalloc(destination.id_low, source.id_low);
  copy_vector_noalloc(destination.position_x_m, source.position_x_m);
  copy_vector_noalloc(destination.position_y_m, source.position_y_m);
  copy_vector_noalloc(destination.position_z_m, source.position_z_m);
  copy_vector_noalloc(destination.velocity_x_m_per_s,
                      source.velocity_x_m_per_s);
  copy_vector_noalloc(destination.velocity_y_m_per_s,
                      source.velocity_y_m_per_s);
  copy_vector_noalloc(destination.velocity_z_m_per_s,
                      source.velocity_z_m_per_s);
  copy_vector_noalloc(destination.droplet_mass_kg,
                      source.droplet_mass_kg);
  copy_vector_noalloc(destination.droplet_diameter_m,
                      source.droplet_diameter_m);
  copy_vector_noalloc(destination.multiplicity, source.multiplicity);
  copy_vector_noalloc(destination.temperature_k, source.temperature_k);
  copy_vector_noalloc(destination.liquid_material_fingerprint,
                      source.liquid_material_fingerprint);
  copy_vector_noalloc(destination.owner_global_cell,
                      source.owner_global_cell);
  copy_vector_noalloc(destination.age_s, source.age_s);
  copy_vector_noalloc(destination.tab_deformation,
                      source.tab_deformation);
  copy_vector_noalloc(destination.tab_deformation_rate_per_s,
                      source.tab_deformation_rate_per_s);
}

SprayParcelState parcel_at(const ParcelSoASnapshot& arrays,
                           std::size_t index) noexcept {
  SprayParcelState parcel;
  parcel.id = id_at(arrays, index);
  parcel.position_m = {arrays.position_x_m[index],
                       arrays.position_y_m[index],
                       arrays.position_z_m[index]};
  parcel.velocity_m_per_s = {arrays.velocity_x_m_per_s[index],
                             arrays.velocity_y_m_per_s[index],
                             arrays.velocity_z_m_per_s[index]};
  parcel.droplet_mass_kg = arrays.droplet_mass_kg[index];
  parcel.droplet_diameter_m = arrays.droplet_diameter_m[index];
  parcel.multiplicity = arrays.multiplicity[index];
  parcel.temperature_k = arrays.temperature_k[index];
  parcel.liquid_material_fingerprint =
      arrays.liquid_material_fingerprint[index];
  parcel.owner_global_cell = arrays.owner_global_cell[index];
  parcel.age_s = arrays.age_s[index];
  return parcel;
}

void write_parcel(ParcelSoASnapshot& arrays, std::size_t index,
                  const SprayParcelState& parcel) noexcept {
  arrays.id_high[index] = parcel.id.high;
  arrays.id_low[index] = parcel.id.low;
  arrays.position_x_m[index] = parcel.position_m[0U];
  arrays.position_y_m[index] = parcel.position_m[1U];
  arrays.position_z_m[index] = parcel.position_m[2U];
  arrays.velocity_x_m_per_s[index] = parcel.velocity_m_per_s[0U];
  arrays.velocity_y_m_per_s[index] = parcel.velocity_m_per_s[1U];
  arrays.velocity_z_m_per_s[index] = parcel.velocity_m_per_s[2U];
  arrays.droplet_mass_kg[index] = parcel.droplet_mass_kg;
  arrays.droplet_diameter_m[index] = parcel.droplet_diameter_m;
  arrays.multiplicity[index] = parcel.multiplicity;
  arrays.temperature_k[index] = parcel.temperature_k;
  arrays.liquid_material_fingerprint[index] =
      parcel.liquid_material_fingerprint;
  arrays.owner_global_cell[index] = parcel.owner_global_cell;
  arrays.age_s[index] = parcel.age_s;
}

void copy_entry(ParcelSoASnapshot& arrays, std::size_t destination,
                std::size_t source) noexcept {
  write_parcel(arrays, destination, parcel_at(arrays, source));
  arrays.tab_deformation[destination] = arrays.tab_deformation[source];
  arrays.tab_deformation_rate_per_s[destination] =
      arrays.tab_deformation_rate_per_s[source];
}

std::size_t lower_bound_index(const ParcelSoASnapshot& arrays,
                              ParcelId id) noexcept {
  std::size_t first = 0U;
  std::size_t count = arrays.size();
  while (count != 0U) {
    const std::size_t step = count / 2U;
    const std::size_t middle = first + step;
    if (id_at(arrays, middle) < id) {
      first = middle + 1U;
      count -= step + 1U;
    } else {
      count = step;
    }
  }
  return first;
}

bool valid_ordered_snapshot(const ParcelSoASnapshot& snapshot) noexcept {
  if (!snapshot.consistent()) {
    return false;
  }
  for (std::size_t index = 0U; index < snapshot.size(); ++index) {
    if (validate_parcel_state(parcel_at(snapshot, index)) !=
            ParcelStateStatus::success ||
        !std::isfinite(snapshot.tab_deformation[index]) ||
        !std::isfinite(snapshot.tab_deformation_rate_per_s[index])) {
      return false;
    }
    if (index != 0U && !(id_at(snapshot, index - 1U) <
                         id_at(snapshot, index))) {
      return false;
    }
  }
  return true;
}

bool valid_injector_spec(const InjectorSpec& spec,
                         InjectorCommittedState initial) noexcept {
  const double axis_norm_squared =
      spec.axis[0U] * spec.axis[0U] + spec.axis[1U] * spec.axis[1U] +
      spec.axis[2U] * spec.axis[2U];
  const bool valid_shape = spec.shape == InjectorShape::point ||
                           spec.shape == InjectorShape::cone;
  return valid_shape && finite_vector(spec.origin_m) &&
         finite_vector(spec.axis) && std::isfinite(axis_norm_squared) &&
         axis_norm_squared > 0.0 &&
         std::isfinite(spec.cone_half_angle_rad) &&
         spec.cone_half_angle_rad >= 0.0 &&
         spec.cone_half_angle_rad <= kPi &&
         std::isfinite(spec.injection_speed_m_per_s) &&
         spec.injection_speed_m_per_s >= 0.0 &&
         std::isfinite(spec.mass_flow_rate_kg_per_s) &&
         spec.mass_flow_rate_kg_per_s >= 0.0 &&
         std::isfinite(spec.represented_mass_per_parcel_kg) &&
         spec.represented_mass_per_parcel_kg > 0.0 &&
         std::isfinite(spec.droplet_mass_kg) &&
         spec.droplet_mass_kg > 0.0 &&
         std::isfinite(spec.droplet_diameter_m) &&
         spec.droplet_diameter_m > 0.0 &&
         std::isfinite(spec.temperature_k) && spec.temperature_k > 0.0 &&
         spec.liquid_material_fingerprint != 0U &&
         std::isfinite(initial.residual_mass_kg) &&
         initial.residual_mass_kg >= 0.0 &&
         initial.residual_mass_kg < spec.represented_mass_per_parcel_kg;
}

Vector3 normalized(Vector3 value) noexcept {
  const double magnitude =
      std::sqrt(value[0U] * value[0U] + value[1U] * value[1U] +
                value[2U] * value[2U]);
  return {value[0U] / magnitude, value[1U] / magnitude,
          value[2U] / magnitude};
}

Vector3 cross(Vector3 left, Vector3 right) noexcept {
  return {left[1U] * right[2U] - left[2U] * right[1U],
          left[2U] * right[0U] - left[0U] * right[2U],
          left[0U] * right[1U] - left[1U] * right[0U]};
}

Vector3 injection_direction(const InjectorSpec& spec,
                            std::uint64_t accepted_step,
                            std::uint64_t ordinal) noexcept {
  const Vector3 axis = normalized(spec.axis);
  if (spec.shape == InjectorShape::point ||
      spec.cone_half_angle_rad == 0.0) {
    return axis;
  }

  const ParcelRandomAddress address{
      spec.seed, accepted_step, spec.injector_id, ordinal,
      ParcelRandomPurpose::injection_direction};
  const double polar_uniform = parcel_uniform_01(address, 0U);
  const double azimuth_uniform = parcel_uniform_01(address, 1U);
  const double minimum_cosine = std::cos(spec.cone_half_angle_rad);
  const double cosine = 1.0 - polar_uniform * (1.0 - minimum_cosine);
  const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
  const double azimuth = 2.0 * kPi * azimuth_uniform;
  const Vector3 reference =
      std::abs(axis[2U]) < 0.9 ? Vector3{0.0, 0.0, 1.0}
                              : Vector3{1.0, 0.0, 0.0};
  const Vector3 tangent = normalized(cross(reference, axis));
  const Vector3 bitangent = cross(axis, tangent);
  return {cosine * axis[0U] +
              sine * (std::cos(azimuth) * tangent[0U] +
                      std::sin(azimuth) * bitangent[0U]),
          cosine * axis[1U] +
              sine * (std::cos(azimuth) * tangent[1U] +
                      std::sin(azimuth) * bitangent[1U]),
          cosine * axis[2U] +
              sine * (std::cos(azimuth) * tangent[2U] +
                      std::sin(azimuth) * bitangent[2U])};
}

InjectionReport base_report(const InjectorSpec& spec,
                            InjectorCommittedState committed,
                            std::uint64_t accepted_step) noexcept {
  InjectionReport report;
  report.accepted_step = accepted_step;
  report.injector_id = spec.injector_id;
  report.first_ordinal = committed.next_ordinal;
  report.next_ordinal = committed.next_ordinal;
  report.residual_mass_before_kg = committed.residual_mass_kg;
  report.residual_mass_after_kg = committed.residual_mass_kg;
  return report;
}

} // namespace

bool ParcelSoASnapshot::consistent() const noexcept {
  const std::size_t count = id_high.size();
  return id_low.size() == count && position_x_m.size() == count &&
         position_y_m.size() == count && position_z_m.size() == count &&
         velocity_x_m_per_s.size() == count &&
         velocity_y_m_per_s.size() == count &&
         velocity_z_m_per_s.size() == count &&
         droplet_mass_kg.size() == count &&
         droplet_diameter_m.size() == count && multiplicity.size() == count &&
         temperature_k.size() == count &&
         liquid_material_fingerprint.size() == count &&
         owner_global_cell.size() == count && age_s.size() == count &&
         tab_deformation.size() == count &&
         tab_deformation_rate_per_s.size() == count;
}

Status ParcelContainer::reserve(std::size_t maximum_parcels) noexcept {
  if (trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_already_active);
  }
  if (maximum_parcels < committed_.size()) {
    return failure(StatusCode::invalid_plan,
                   ParcelOperationDetail::capacity_exceeded);
  }
  try {
    reserve_arrays(committed_, maximum_parcels);
    reserve_arrays(trial_, maximum_parcels);
  } catch (const std::bad_alloc&) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::capacity_exceeded);
  } catch (...) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::capacity_exceeded);
  }
  capacity_limit_ = maximum_parcels;
  return {};
}

bool ParcelContainer::committed_at(std::size_t index,
                                   SprayParcelState& out) const noexcept {
  if (index >= committed_.size()) {
    return false;
  }
  out = parcel_at(committed_, index);
  return true;
}

bool ParcelContainer::trial_at(std::size_t index,
                               SprayParcelState& out) const noexcept {
  if (!trial_active_ || index >= trial_.size()) {
    return false;
  }
  out = parcel_at(trial_, index);
  return true;
}

bool ParcelContainer::committed_tab_state(
    ParcelId id, double& deformation,
    double& deformation_rate_per_s) const noexcept {
  const std::size_t index = lower_bound_index(committed_, id);
  if (index == committed_.size() || id_at(committed_, index) != id) {
    return false;
  }
  deformation = committed_.tab_deformation[index];
  deformation_rate_per_s = committed_.tab_deformation_rate_per_s[index];
  return true;
}

bool ParcelContainer::trial_tab_state(
    ParcelId id, double& deformation,
    double& deformation_rate_per_s) const noexcept {
  if (!trial_active_) {
    return false;
  }
  const std::size_t index = lower_bound_index(trial_, id);
  if (index == trial_.size() || id_at(trial_, index) != id) {
    return false;
  }
  deformation = trial_.tab_deformation[index];
  deformation_rate_per_s = trial_.tab_deformation_rate_per_s[index];
  return true;
}

Status ParcelContainer::begin_trial() noexcept {
  if (trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_already_active);
  }
  if (!arrays_have_capacity(trial_, committed_.size())) {
    return failure(StatusCode::invalid_plan,
                   ParcelOperationDetail::capacity_exceeded);
  }
  copy_arrays_noalloc(trial_, committed_);
  trial_failure_ = {};
  trial_active_ = true;
  return {};
}

Status ParcelContainer::latch(Status failure_status) noexcept {
  if (trial_failure_) {
    trial_failure_ = failure_status;
  }
  return trial_failure_;
}

Status ParcelContainer::stage_add(
    const SprayParcelState& parcel) noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  if (!trial_failure_) {
    return trial_failure_;
  }
  if (validate_parcel_state(parcel) != ParcelStateStatus::success) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::invalid_parcel));
  }
  const std::size_t insertion = lower_bound_index(trial_, parcel.id);
  if (insertion != trial_.size() && id_at(trial_, insertion) == parcel.id) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::duplicate_id));
  }
  if (trial_.size() >= capacity_limit_ ||
      !arrays_have_capacity(trial_, trial_.size() + 1U)) {
    return latch(failure(StatusCode::invalid_plan,
                         ParcelOperationDetail::capacity_exceeded));
  }

  const std::size_t old_size = trial_.size();
  resize_arrays(trial_, old_size + 1U);
  for (std::size_t index = old_size; index > insertion; --index) {
    copy_entry(trial_, index, index - 1U);
  }
  write_parcel(trial_, insertion, parcel);
  trial_.tab_deformation[insertion] = 0.0;
  trial_.tab_deformation_rate_per_s[insertion] = 0.0;
  return {};
}

Status ParcelContainer::stage_update(
    const SprayParcelState& parcel) noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  if (!trial_failure_) {
    return trial_failure_;
  }
  if (validate_parcel_state(parcel) != ParcelStateStatus::success) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::invalid_parcel));
  }
  const std::size_t index = lower_bound_index(trial_, parcel.id);
  if (index == trial_.size() || id_at(trial_, index) != parcel.id) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::parcel_not_found));
  }
  write_parcel(trial_, index, parcel);
  return {};
}

Status ParcelContainer::stage_tab_state(
    ParcelId id, double deformation,
    double deformation_rate_per_s) noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  if (!trial_failure_) {
    return trial_failure_;
  }
  if (!valid_id(id) || !std::isfinite(deformation) ||
      !std::isfinite(deformation_rate_per_s)) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::invalid_tab_state));
  }
  const std::size_t index = lower_bound_index(trial_, id);
  if (index == trial_.size() || id_at(trial_, index) != id) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::parcel_not_found));
  }
  trial_.tab_deformation[index] = deformation;
  trial_.tab_deformation_rate_per_s[index] = deformation_rate_per_s;
  return {};
}

Status ParcelContainer::stage_remove(ParcelId id) noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  if (!trial_failure_) {
    return trial_failure_;
  }
  if (!valid_id(id)) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::invalid_parcel));
  }
  const std::size_t removal = lower_bound_index(trial_, id);
  if (removal == trial_.size() || id_at(trial_, removal) != id) {
    return latch(failure(StatusCode::invalid_case,
                         ParcelOperationDetail::parcel_not_found));
  }
  for (std::size_t index = removal + 1U; index < trial_.size(); ++index) {
    copy_entry(trial_, index - 1U, index);
  }
  resize_arrays(trial_, trial_.size() - 1U);
  return {};
}

Status ParcelContainer::preflight_commit() const noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  return trial_failure_;
}

void ParcelContainer::publish_preflighted_trial() noexcept {
  assert(trial_active_ && static_cast<bool>(trial_failure_));
  using std::swap;
  swap(committed_, trial_);
  resize_arrays(trial_, 0U);
  trial_active_ = false;
  trial_failure_ = {};
}

Status ParcelContainer::commit_trial() noexcept {
  const Status ready = preflight_commit();
  if (!ready) {
    if (trial_active_) {
      resize_arrays(trial_, 0U);
      trial_active_ = false;
      trial_failure_ = {};
    }
    return ready;
  }
  publish_preflighted_trial();
  return {};
}

Status ParcelContainer::rollback_trial() noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  resize_arrays(trial_, 0U);
  trial_active_ = false;
  trial_failure_ = {};
  return {};
}

Status ParcelContainer::snapshot_committed(
    ParcelSoASnapshot& out) const noexcept {
  try {
    ParcelSoASnapshot candidate = committed_;
    out = std::move(candidate);
  } catch (const std::bad_alloc&) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::invalid_snapshot);
  } catch (...) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::invalid_snapshot);
  }
  return {};
}

Status ParcelContainer::restore_committed(
    const ParcelSoASnapshot& snapshot) noexcept {
  if (trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_already_active);
  }
  if (!valid_ordered_snapshot(snapshot)) {
    return failure(StatusCode::invalid_case,
                   ParcelOperationDetail::invalid_snapshot);
  }
  if (snapshot.size() > capacity_limit_ ||
      !arrays_have_capacity(committed_, snapshot.size())) {
    return failure(StatusCode::invalid_plan,
                   ParcelOperationDetail::capacity_exceeded);
  }
  copy_arrays_noalloc(committed_, snapshot);
  resize_arrays(trial_, 0U);
  trial_failure_ = {};
  return {};
}

Status DeterministicInjector::reserve(
    std::size_t maximum_candidates) noexcept {
  if (trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_already_active);
  }
  try {
    reserve_arrays(candidates_, maximum_candidates);
    candidate_id_scratch_.reserve(maximum_candidates);
  } catch (const std::bad_alloc&) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::capacity_exceeded);
  } catch (...) {
    return failure(StatusCode::allocation_failure,
                   ParcelOperationDetail::capacity_exceeded);
  }
  candidate_capacity_ = maximum_candidates;
  return {};
}

Status DeterministicInjector::configure(
    const InjectorSpec& spec,
    InjectorCommittedState initial_state) noexcept {
  if (trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_already_active);
  }
  if (!valid_injector_spec(spec, initial_state)) {
    return failure(StatusCode::invalid_case,
                   ParcelOperationDetail::invalid_injector);
  }
  spec_ = spec;
  committed_ = initial_state;
  trial_state_ = initial_state;
  resize_arrays(candidates_, 0U);
  candidate_id_scratch_.clear();
  configured_ = true;
  report_ = base_report(spec_, committed_, 0U);
  return {};
}

InjectionReport DeterministicInjector::begin_trial(
    std::uint64_t accepted_step, double duration_s) noexcept {
  if (!configured_) {
    InjectionReport report;
    report.status = InjectionStatus::not_configured;
    report.accepted_step = accepted_step;
    return report;
  }
  if (trial_active_) {
    InjectionReport report = base_report(spec_, committed_, accepted_step);
    report.status = InjectionStatus::trial_already_active;
    return report;
  }

  resize_arrays(candidates_, 0U);
  candidate_id_scratch_.clear();
  report_ = base_report(spec_, committed_, accepted_step);
  if (!std::isfinite(duration_s) || duration_s < 0.0) {
    report_.status = InjectionStatus::invalid_input;
    return report_;
  }

  if (duration_s == 0.0 || spec_.mass_flow_rate_kg_per_s == 0.0) {
    report_.status = InjectionStatus::success;
    trial_state_ = committed_;
    trial_active_ = true;
    return report_;
  }

  const double requested_mass = spec_.mass_flow_rate_kg_per_s * duration_s;
  const double available_mass = committed_.residual_mass_kg + requested_mass;
  if (!std::isfinite(requested_mass) || !std::isfinite(available_mass)) {
    report_.status = InjectionStatus::non_finite_output;
    return report_;
  }
  report_.requested_mass_kg = requested_mass;
  const double count_as_double =
      std::floor(available_mass / spec_.represented_mass_per_parcel_kg);
  if (!std::isfinite(count_as_double) || count_as_double < 0.0 ||
      count_as_double > kMaximumExactlyRepresentableCount ||
      count_as_double > static_cast<double>(candidate_capacity_)) {
    report_.status = InjectionStatus::capacity_exceeded;
    return report_;
  }
  const std::size_t count = static_cast<std::size_t>(count_as_double);
  if (count > candidate_capacity_ ||
      !arrays_have_capacity(candidates_, count)) {
    report_.status = InjectionStatus::capacity_exceeded;
    return report_;
  }
  if (count >
      std::numeric_limits<std::uint64_t>::max() - committed_.next_ordinal) {
    report_.status = InjectionStatus::counter_overflow;
    return report_;
  }

  const std::uint64_t count_u64 = static_cast<std::uint64_t>(count);
  const double injected_mass =
      static_cast<double>(count_u64) *
      spec_.represented_mass_per_parcel_kg;
  double residual_mass = available_mass - injected_mass;
  const double mass_scale =
      std::max({std::abs(available_mass), std::abs(injected_mass),
                std::abs(spec_.represented_mass_per_parcel_kg)});
  const double residual_tolerance =
      8.0 * std::numeric_limits<double>::epsilon() * mass_scale;
  if (residual_mass < 0.0 && residual_mass >= -residual_tolerance) {
    residual_mass = 0.0;
  }
  if (!std::isfinite(injected_mass) || !std::isfinite(residual_mass) ||
      residual_mass < 0.0 ||
      residual_mass >= spec_.represented_mass_per_parcel_kg) {
    report_.status = InjectionStatus::non_finite_output;
    return report_;
  }

  const double multiplicity =
      spec_.represented_mass_per_parcel_kg / spec_.droplet_mass_kg;
  if (!std::isfinite(multiplicity) || multiplicity <= 0.0) {
    report_.status = InjectionStatus::non_finite_output;
    return report_;
  }

  resize_arrays(candidates_, count);
  candidate_id_scratch_.resize(count);
  for (std::size_t index = 0U; index < count; ++index) {
    const std::uint64_t ordinal =
        committed_.next_ordinal + static_cast<std::uint64_t>(index);
    const ParcelRandomAddress id_address{
        spec_.seed, accepted_step, spec_.injector_id, ordinal,
        ParcelRandomPurpose::stable_id};
    const Vector3 direction =
        injection_direction(spec_, accepted_step, ordinal);
    SprayParcelState parcel;
    parcel.id = make_stable_parcel_id(id_address);
    parcel.position_m = spec_.origin_m;
    parcel.velocity_m_per_s = {
        direction[0U] * spec_.injection_speed_m_per_s,
        direction[1U] * spec_.injection_speed_m_per_s,
        direction[2U] * spec_.injection_speed_m_per_s};
    parcel.droplet_mass_kg = spec_.droplet_mass_kg;
    parcel.droplet_diameter_m = spec_.droplet_diameter_m;
    parcel.multiplicity = multiplicity;
    parcel.temperature_k = spec_.temperature_k;
    parcel.liquid_material_fingerprint =
        spec_.liquid_material_fingerprint;
    parcel.owner_global_cell = spec_.owner_global_cell;
    parcel.age_s = 0.0;
    if (validate_parcel_state(parcel) != ParcelStateStatus::success) {
      resize_arrays(candidates_, 0U);
      candidate_id_scratch_.clear();
      report_.status = InjectionStatus::non_finite_output;
      return report_;
    }
    write_parcel(candidates_, index, parcel);
    candidate_id_scratch_[index] = parcel.id;
  }
  std::sort(candidate_id_scratch_.begin(), candidate_id_scratch_.end());
  for (std::size_t index = 1U; index < candidate_id_scratch_.size();
       ++index) {
    if (candidate_id_scratch_[index] == candidate_id_scratch_[index - 1U]) {
      resize_arrays(candidates_, 0U);
      candidate_id_scratch_.clear();
      report_.status = InjectionStatus::non_finite_output;
      return report_;
    }
  }

  trial_state_ = {residual_mass,
                  committed_.next_ordinal + count_u64};
  report_.status = InjectionStatus::success;
  report_.next_ordinal = trial_state_.next_ordinal;
  report_.parcel_count = count;
  report_.injected_mass_kg = injected_mass;
  report_.residual_mass_after_kg = residual_mass;
  trial_active_ = true;
  return report_;
}

Status DeterministicInjector::preflight_commit() const noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  return {};
}

void DeterministicInjector::publish_preflighted_trial() noexcept {
  assert(trial_active_);
  committed_ = trial_state_;
  resize_arrays(candidates_, 0U);
  candidate_id_scratch_.clear();
  trial_active_ = false;
}

Status DeterministicInjector::commit_trial() noexcept {
  const Status ready = preflight_commit();
  if (!ready) {
    return ready;
  }
  publish_preflighted_trial();
  return {};
}

Status DeterministicInjector::rollback_trial() noexcept {
  if (!trial_active_) {
    return failure(StatusCode::rejected_step,
                   ParcelOperationDetail::trial_not_active);
  }
  resize_arrays(candidates_, 0U);
  candidate_id_scratch_.clear();
  trial_state_ = committed_;
  trial_active_ = false;
  return {};
}

bool DeterministicInjector::candidate_at(
    std::size_t index, SprayParcelState& out) const noexcept {
  if (!trial_active_ || index >= candidates_.size()) {
    return false;
  }
  out = parcel_at(candidates_, index);
  return true;
}

} // namespace hundun::v04::spray::detail

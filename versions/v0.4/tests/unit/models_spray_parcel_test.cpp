// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_parcel_detail.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

std::size_t hot_allocation_count = 0U;
bool count_hot_allocations = false;

} // namespace

void* operator new(std::size_t size) {
  if (count_hot_allocations) {
    ++hot_allocation_count;
  }
  if (void* storage = std::malloc(size == 0U ? 1U : size)) {
    return storage;
  }
  throw std::bad_alloc{};
}

void* operator new[](std::size_t size) {
  return ::operator new(size);
}

void operator delete(void* storage) noexcept { std::free(storage); }
void operator delete[](void* storage) noexcept { std::free(storage); }
void operator delete(void* storage, std::size_t) noexcept {
  std::free(storage);
}
void operator delete[](void* storage, std::size_t) noexcept {
  std::free(storage);
}

namespace {

using hundun::v04::Status;
using hundun::v04::StatusCode;
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

SprayParcelState parcel(std::uint64_t high, std::uint64_t low,
                        double marker) {
  SprayParcelState value;
  value.id = {high, low};
  value.position_m = {marker + 0.1, marker + 0.2, marker + 0.3};
  value.velocity_m_per_s = {marker + 1.1, marker + 1.2, marker + 1.3};
  value.droplet_mass_kg = marker + 2.1;
  value.droplet_diameter_m = marker + 2.2;
  value.multiplicity = marker + 2.3;
  value.temperature_k = marker + 300.0;
  value.liquid_material_fingerprint = high + 100U;
  value.owner_global_cell = low + 200U;
  value.age_s = marker + 0.4;
  return value;
}

bool same_parcel(const SprayParcelState& left,
                 const SprayParcelState& right) {
  return left.id == right.id && left.position_m == right.position_m &&
         left.velocity_m_per_s == right.velocity_m_per_s &&
         left.droplet_mass_kg == right.droplet_mass_kg &&
         left.droplet_diameter_m == right.droplet_diameter_m &&
         left.multiplicity == right.multiplicity &&
         left.temperature_k == right.temperature_k &&
         left.liquid_material_fingerprint ==
             right.liquid_material_fingerprint &&
         left.owner_global_cell == right.owner_global_cell &&
         left.age_s == right.age_s;
}

template <class T>
bool same_vector_bytes(const std::vector<T>& left,
                       const std::vector<T>& right) {
  static_assert(std::is_trivially_copyable<T>::value,
                "parcel snapshots contain primitive arrays");
  return left.size() == right.size() &&
         (left.empty() ||
          std::memcmp(left.data(), right.data(), left.size() * sizeof(T)) ==
              0);
}

bool same_snapshot_bytes(const ParcelSoASnapshot& left,
                         const ParcelSoASnapshot& right) {
  return same_vector_bytes(left.id_high, right.id_high) &&
         same_vector_bytes(left.id_low, right.id_low) &&
         same_vector_bytes(left.position_x_m, right.position_x_m) &&
         same_vector_bytes(left.position_y_m, right.position_y_m) &&
         same_vector_bytes(left.position_z_m, right.position_z_m) &&
         same_vector_bytes(left.velocity_x_m_per_s,
                           right.velocity_x_m_per_s) &&
         same_vector_bytes(left.velocity_y_m_per_s,
                           right.velocity_y_m_per_s) &&
         same_vector_bytes(left.velocity_z_m_per_s,
                           right.velocity_z_m_per_s) &&
         same_vector_bytes(left.droplet_mass_kg, right.droplet_mass_kg) &&
         same_vector_bytes(left.droplet_diameter_m,
                           right.droplet_diameter_m) &&
         same_vector_bytes(left.multiplicity, right.multiplicity) &&
         same_vector_bytes(left.temperature_k, right.temperature_k) &&
         same_vector_bytes(left.liquid_material_fingerprint,
                           right.liquid_material_fingerprint) &&
         same_vector_bytes(left.owner_global_cell,
                           right.owner_global_cell) &&
         same_vector_bytes(left.age_s, right.age_s) &&
         same_vector_bytes(left.tab_deformation, right.tab_deformation) &&
         same_vector_bytes(left.tab_deformation_rate_per_s,
                           right.tab_deformation_rate_per_s);
}

bool snapshot(const ParcelContainer& container, ParcelSoASnapshot& out) {
  return static_cast<bool>(container.snapshot_committed(out));
}

bool test_container_transaction_and_snapshot() {
  static_assert(!std::is_move_constructible<ParcelContainer>::value,
                "moving a transactional container would invalidate its layers");
  ParcelContainer container;
  bool passed = true;
  passed &= expect(static_cast<bool>(container.reserve(4U)),
                   "container reserves both state layers");

  hot_allocation_count = 0U;
  count_hot_allocations = true;
  const Status begin = container.begin_trial();
  const Status add_3 = container.stage_add(parcel(2U, 9U, 3.0));
  const Status add_1 = container.stage_add(parcel(1U, 9U, 1.0));
  const Status add_2 = container.stage_add(parcel(2U, 8U, 2.0));
  const Status tab_3 = container.stage_tab_state({2U, 9U}, 0.3, 30.0);
  const Status tab_1 = container.stage_tab_state({1U, 9U}, 0.1, 10.0);
  const Status tab_2 = container.stage_tab_state({2U, 8U}, 0.2, 20.0);
  const Status ready = container.preflight_commit();
  const Status ready_again = container.preflight_commit();
  if (ready && ready_again) {
    container.publish_preflighted_trial();
  }
  count_hot_allocations = false;
  passed &= expect(static_cast<bool>(begin) && static_cast<bool>(add_3) &&
                       static_cast<bool>(add_1) &&
                       static_cast<bool>(add_2) &&
                       static_cast<bool>(tab_3) &&
                       static_cast<bool>(tab_1) &&
                       static_cast<bool>(tab_2) &&
                       static_cast<bool>(ready) &&
                       static_cast<bool>(ready_again) &&
                       !container.trial_active() &&
                       hot_allocation_count == 0U,
                   "read-only preflight and no-fail publish allocate nothing");
  passed &= expect(container.committed_size() == 3U,
                   "three parcels are committed");

  SprayParcelState observed;
  passed &= expect(container.committed_at(0U, observed) &&
                       observed.id == ParcelId{1U, 9U} &&
                       container.committed_at(1U, observed) &&
                       observed.id == ParcelId{2U, 8U} &&
                       container.committed_at(2U, observed) &&
                       observed.id == ParcelId{2U, 9U},
                   "committed SoA is ordered by the complete 128-bit ID");
  double deformation = 0.0;
  double deformation_rate = 0.0;
  passed &= expect(container.committed_tab_state(
                       {1U, 9U}, deformation, deformation_rate) &&
                       deformation == 0.1 && deformation_rate == 10.0 &&
                       container.committed_tab_state(
                           {2U, 8U}, deformation, deformation_rate) &&
                       deformation == 0.2 && deformation_rate == 20.0,
                   "TAB continuation state follows its stable parcel ID");

  ParcelSoASnapshot baseline;
  passed &= expect(snapshot(container, baseline),
                   "committed snapshot is available");

  passed &= expect(static_cast<bool>(container.begin_trial()),
                   "duplicate attempt begins");
  SprayParcelState changed = parcel(2U, 8U, 22.0);
  passed &= expect(static_cast<bool>(container.stage_update(changed)),
                   "an update can be staged before a later failure");
  passed &= expect(container.trial_tab_state(
                       {2U, 8U}, deformation, deformation_rate) &&
                       deformation == 0.2 && deformation_rate == 20.0,
                   "updating public parcel values preserves TAB state");
  const Status duplicate = container.stage_add(parcel(2U, 8U, 2.0));
  const Status failed_commit = container.commit_trial();
  ParcelSoASnapshot after_duplicate;
  passed &= expect(duplicate.code == StatusCode::invalid_case &&
                       duplicate.detail == static_cast<std::uint32_t>(
                                               ParcelOperationDetail::duplicate_id) &&
                       failed_commit.code == duplicate.code &&
                       failed_commit.detail == duplicate.detail &&
                       snapshot(container, after_duplicate) &&
                       same_snapshot_bytes(after_duplicate, baseline),
                   "duplicate ID poisons the trial and publishes nothing");

  passed &= expect(static_cast<bool>(container.begin_trial()) &&
                       static_cast<bool>(container.stage_remove({1U, 9U})) &&
                       static_cast<bool>(container.stage_update(changed)) &&
                       static_cast<bool>(container.stage_tab_state(
                           {2U, 8U}, 8.0, 80.0)) &&
                       static_cast<bool>(container.rollback_trial()),
                   "remove/update trial can be rolled back");
  ParcelSoASnapshot after_rollback;
  passed &= expect(snapshot(container, after_rollback) &&
                       same_snapshot_bytes(after_rollback, baseline),
                   "rollback preserves every committed parcel byte");

  passed &= expect(static_cast<bool>(container.begin_trial()) &&
                       static_cast<bool>(
                           container.stage_add(parcel(3U, 6U, 4.0))),
                   "fourth parcel fills reserved capacity");
  const Status overflow = container.stage_add(parcel(4U, 5U, 5.0));
  const Status overflow_commit = container.commit_trial();
  ParcelSoASnapshot after_overflow;
  passed &= expect(overflow.code == StatusCode::invalid_plan &&
                       overflow.detail ==
                           static_cast<std::uint32_t>(
                               ParcelOperationDetail::capacity_exceeded) &&
                       overflow_commit.code == overflow.code &&
                       overflow_commit.detail == overflow.detail &&
                       snapshot(container, after_overflow) &&
                       same_snapshot_bytes(after_overflow, baseline),
                   "capacity overflow cannot publish a partial trial");

  ParcelContainer restored;
  ParcelSoASnapshot restored_snapshot;
  passed &= expect(static_cast<bool>(restored.reserve(4U)) &&
                       static_cast<bool>(restored.restore_committed(baseline)) &&
                       snapshot(restored, restored_snapshot) &&
                       same_snapshot_bytes(restored_snapshot, baseline) &&
                       restored.committed_tab_state(
                           {2U, 9U}, deformation, deformation_rate) &&
                       deformation == 0.3 && deformation_rate == 30.0,
                   "SoA snapshot round-trip preserves parcel and TAB state");

  passed &= expect(static_cast<bool>(container.begin_trial()) &&
                       static_cast<bool>(container.stage_remove({1U, 9U})) &&
                       static_cast<bool>(container.stage_update(changed)) &&
                       static_cast<bool>(container.stage_tab_state(
                           {2U, 8U}, 8.0, 80.0)) &&
                       static_cast<bool>(
                           container.stage_add(parcel(3U, 6U, 4.0))) &&
                       static_cast<bool>(container.preflight_commit()),
                   "staged add/update/remove commit together");
  container.publish_preflighted_trial();
  passed &= expect(container.committed_size() == 3U &&
                       container.committed_at(0U, observed) &&
                       same_parcel(observed, changed) &&
                       container.committed_tab_state(
                           {2U, 8U}, deformation, deformation_rate) &&
                       deformation == 8.0 && deformation_rate == 80.0 &&
                       container.committed_tab_state(
                           {2U, 9U}, deformation, deformation_rate) &&
                       deformation == 0.3 && deformation_rate == 30.0 &&
                       container.committed_tab_state(
                           {3U, 6U}, deformation, deformation_rate) &&
                       deformation == 0.0 && deformation_rate == 0.0,
                   "commit keeps shifted TAB state and zeros a new parcel");
  return passed;
}

InjectorSpec base_injector() {
  InjectorSpec spec;
  spec.seed = 777U;
  spec.injector_id = 31U;
  spec.shape = InjectorShape::point;
  spec.origin_m = {1.0, 2.0, 3.0};
  spec.axis = {2.0, 0.0, 0.0};
  spec.cone_half_angle_rad = 0.0;
  spec.injection_speed_m_per_s = 10.0;
  spec.mass_flow_rate_kg_per_s = 2.5;
  spec.represented_mass_per_parcel_kg = 1.0;
  spec.droplet_mass_kg = 0.25;
  spec.droplet_diameter_m = 1.0e-3;
  spec.temperature_k = 345.0;
  spec.liquid_material_fingerprint = 4242U;
  spec.owner_global_cell = 99U;
  return spec;
}

bool same_report(const InjectionReport& left, const InjectionReport& right) {
  return left.status == right.status && left.model_id == right.model_id &&
         left.accepted_step == right.accepted_step &&
         left.injector_id == right.injector_id &&
         left.first_ordinal == right.first_ordinal &&
         left.next_ordinal == right.next_ordinal &&
         left.parcel_count == right.parcel_count &&
         left.requested_mass_kg == right.requested_mass_kg &&
         left.injected_mass_kg == right.injected_mass_kg &&
         left.residual_mass_before_kg == right.residual_mass_before_kg &&
         left.residual_mass_after_kg == right.residual_mass_after_kg;
}

bool test_injector_residual_retry_and_identity() {
  static_assert(!std::is_move_constructible<DeterministicInjector>::value,
                "moving an active injector would invalidate its trial state");
  DeterministicInjector injector;
  const InjectorSpec spec = base_injector();
  bool passed = true;
  passed &= expect(static_cast<bool>(injector.reserve(4U)) &&
                       static_cast<bool>(injector.configure(
                           spec, InjectorCommittedState{0.25, 100U})),
                   "injector configuration is explicit and preallocated");

  hot_allocation_count = 0U;
  count_hot_allocations = true;
  const InjectionReport first = injector.begin_trial(17U, 1.0);
  count_hot_allocations = false;
  passed &= expect(first.succeeded() && first.parcel_count == 2U &&
                       first.requested_mass_kg == 2.5 &&
                       first.injected_mass_kg == 2.0 &&
                       first.residual_mass_before_kg == 0.25 &&
                       first.residual_mass_after_kg == 0.75 &&
                       first.first_ordinal == 100U &&
                       first.next_ordinal == 102U &&
                       hot_allocation_count == 0U,
                   "fractional represented mass is retained without hot allocation");
  passed &= expect(std::abs((first.residual_mass_before_kg +
                             first.requested_mass_kg) -
                            (first.injected_mass_kg +
                             first.residual_mass_after_kg)) <=
                       8.0 * std::numeric_limits<double>::epsilon(),
                   "injected mass plus residual conserves requested mass");
  passed &= expect(injector.committed_state() ==
                       InjectorCommittedState{0.25, 100U},
                   "trial generation does not consume committed state");

  std::vector<SprayParcelState> first_candidates;
  for (std::size_t index = 0U; index < injector.candidate_count(); ++index) {
    SprayParcelState value;
    passed &= expect(injector.candidate_at(index, value),
                     "trial candidate is readable");
    first_candidates.push_back(value);
  }
  passed &= expect(first_candidates.size() == 2U &&
                       first_candidates[0U].id == make_stable_parcel_id(
                           {spec.seed, 17U, spec.injector_id, 100U,
                            ParcelRandomPurpose::stable_id}) &&
                       first_candidates[1U].id == make_stable_parcel_id(
                           {spec.seed, 17U, spec.injector_id, 101U,
                            ParcelRandomPurpose::stable_id}) &&
                       first_candidates[0U].velocity_m_per_s ==
                           Vector3{10.0, 0.0, 0.0} &&
                       first_candidates[0U].multiplicity == 4.0 &&
                       first_candidates[0U].liquid_material_fingerprint ==
                           4242U,
                   "point candidates use stable IDs and configured liquid data");

  passed &= expect(static_cast<bool>(injector.rollback_trial()) &&
                       injector.committed_state() ==
                           InjectorCommittedState{0.25, 100U},
                   "rollback consumes neither residual nor ordinal");
  const InjectionReport retry = injector.begin_trial(17U, 1.0);
  passed &= expect(same_report(first, retry) &&
                       injector.candidate_count() == first_candidates.size(),
                   "retry report is deterministic");
  for (std::size_t index = 0U; index < first_candidates.size(); ++index) {
    SprayParcelState value;
    passed &= expect(injector.candidate_at(index, value) &&
                         same_parcel(value, first_candidates[index]),
                     "retry parcel candidate is bitwise deterministic by value");
  }
  passed &= expect(static_cast<bool>(injector.preflight_commit()),
                   "injector preflight is read-only");
  injector.publish_preflighted_trial();
  passed &= expect(!injector.trial_active() &&
                       injector.committed_state() ==
                           InjectorCommittedState{0.75, 102U},
                   "commit advances residual and ordinal once");
  return passed;
}

bool test_micro_mass_residual_scale() {
  InjectorSpec micro = base_injector();
  micro.mass_flow_rate_kg_per_s = 2.5e-30;
  micro.represented_mass_per_parcel_kg = 1.0e-30;
  micro.droplet_mass_kg = 2.5e-31;
  micro.droplet_diameter_m = 1.0e-9;

  DeterministicInjector injector;
  bool passed = true;
  passed &= expect(static_cast<bool>(injector.reserve(3U)) &&
                       static_cast<bool>(injector.configure(micro)),
                   "microscopic-mass injector configures");
  const InjectionReport report = injector.begin_trial(6U, 1.0);
  const double balance =
      std::abs((report.residual_mass_before_kg + report.requested_mass_kg) -
               (report.injected_mass_kg + report.residual_mass_after_kg));
  const double scale = std::abs(report.requested_mass_kg) +
                       std::abs(report.injected_mass_kg) +
                       std::abs(report.residual_mass_after_kg);
  passed &= expect(report.succeeded() && report.parcel_count == 2U &&
                       report.residual_mass_after_kg > 0.0 &&
                       report.residual_mass_after_kg <
                           micro.represented_mass_per_parcel_kg &&
                       balance <=
                           8.0 * std::numeric_limits<double>::epsilon() * scale,
                   "roundoff budget scales with physical parcel mass");
  passed &= expect(static_cast<bool>(injector.rollback_trial()),
                   "micro-mass trial remains rollback-safe");
  return passed;
}

bool test_cone_zero_flow_and_capacity_failure() {
  bool passed = true;
  InjectorSpec cone = base_injector();
  cone.shape = InjectorShape::cone;
  cone.axis = {0.0, 0.0, 4.0};
  cone.cone_half_angle_rad = 0.25;
  cone.mass_flow_rate_kg_per_s = 2.0;

  DeterministicInjector cone_injector;
  passed &= expect(static_cast<bool>(cone_injector.reserve(2U)) &&
                       static_cast<bool>(cone_injector.configure(cone)),
                   "cone injector configures");
  const InjectionReport cone_report = cone_injector.begin_trial(3U, 1.0);
  passed &= expect(cone_report.succeeded() &&
                       cone_report.parcel_count == 2U,
                   "cone produces deterministic candidates");
  for (std::size_t index = 0U; index < cone_injector.candidate_count();
       ++index) {
    SprayParcelState value;
    passed &= expect(cone_injector.candidate_at(index, value),
                     "cone candidate is readable");
    const double speed = std::sqrt(
        value.velocity_m_per_s[0U] * value.velocity_m_per_s[0U] +
        value.velocity_m_per_s[1U] * value.velocity_m_per_s[1U] +
        value.velocity_m_per_s[2U] * value.velocity_m_per_s[2U]);
    const double cosine = value.velocity_m_per_s[2U] / speed;
    passed &= expect(std::abs(speed - cone.injection_speed_m_per_s) <=
                             16.0 * std::numeric_limits<double>::epsilon() *
                                 cone.injection_speed_m_per_s &&
                         cosine >= std::cos(cone.cone_half_angle_rad) -
                                       16.0 *
                                           std::numeric_limits<double>::epsilon() &&
                         cosine <= 1.0 +
                                       16.0 *
                                           std::numeric_limits<double>::epsilon(),
                     "cone direction remains inside its angular bound");
  }
  passed &= expect(static_cast<bool>(cone_injector.rollback_trial()),
                   "cone candidate can be discarded");

  InjectorSpec stopped = base_injector();
  stopped.mass_flow_rate_kg_per_s = 0.0;
  DeterministicInjector zero_flow;
  passed &= expect(static_cast<bool>(zero_flow.reserve(2U)) &&
                       static_cast<bool>(zero_flow.configure(
                           stopped, InjectorCommittedState{0.25, 9U})),
                   "zero-flow injector configures");
  const InjectionReport zero = zero_flow.begin_trial(4U, 7.0);
  passed &= expect(zero.succeeded() && zero.parcel_count == 0U &&
                       zero.requested_mass_kg == 0.0 &&
                       zero.injected_mass_kg == 0.0 &&
                       zero.residual_mass_before_kg == 0.25 &&
                       zero.residual_mass_after_kg == 0.25 &&
                       static_cast<bool>(zero_flow.commit_trial()) &&
                       zero_flow.committed_state() ==
                           InjectorCommittedState{0.25, 9U},
                   "no mass flow is an exact zero and preserves residual state");

  InjectorSpec too_many = base_injector();
  too_many.mass_flow_rate_kg_per_s = 3.0;
  DeterministicInjector bounded;
  passed &= expect(static_cast<bool>(bounded.reserve(1U)) &&
                       static_cast<bool>(bounded.configure(
                           too_many, InjectorCommittedState{0.0, 40U})),
                   "bounded injector configures");
  const InjectionReport overflow = bounded.begin_trial(5U, 1.0);
  passed &= expect(overflow.status == InjectionStatus::capacity_exceeded &&
                       !bounded.trial_active() &&
                       bounded.candidate_count() == 0U &&
                       bounded.committed_state() ==
                           InjectorCommittedState{0.0, 40U},
                   "candidate overflow is explicit and cannot consume state");
  return passed;
}

} // namespace

int main() {
  bool passed = test_container_transaction_and_snapshot();
  passed &= test_injector_residual_retry_and_identity();
  passed &= test_micro_mass_residual_scale();
  passed &= test_cone_zero_flow_and_capacity_failure();
  return passed ? 0 : 1;
}

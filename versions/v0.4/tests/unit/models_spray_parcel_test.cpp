// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "core_spray_history_detail.hpp"
#include "core_tcr_history_detail.hpp"
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

bool test_native_spray_history(bool coupled_tcr) {
  using namespace hundun::v04;
  hundun::v04::detail::ProductSprayHistory history;
  hundun::v04::detail::ProductTcrHistory tcr;
  if (coupled_tcr)
    tcr.configure(9001, 8, -1);
  const auto stage_tcr = [&] {
    if (!coupled_tcr)
      return true;
    for (std::size_t cell = 0; cell < 8; ++cell) {
      const auto &old = tcr.accepted(cell);
      tcr::detail::TrialRequest request;
      request.expected_revision = old.revision;
      request.mode = tcr::detail::Mode::experimental;
      request.initialization_sign = -1;
      request.mapping = {tcr::detail::Status::success,
                         {.25, 1},
                         tcr::detail::kReactantMoleFractionMappingIdentity};
      const auto trial = tcr::detail::prepare(old, request);
      if (!tcr.stage(cell, trial, 0))
        return false;
    }
    tcr.seal();
    return true;
  };
  auto spec = base_injector();
  spec.owner_global_cell = 7;
  const InjectorCommittedState original{.25, UINT64_C(9007199254740993)};
  DeterministicInjector injector;
  if (!injector.reserve(4) || !injector.configure(spec, original))
    return false;
  DeterministicInjector *injections[]{&injector};
  const Int3 global{2, 2, 2};
  const MeshPatch patch{{0, 0, 0}, global, {1, 1, 1}, {0, 0, 0}};
  if (!history.configure(8001, patch, global, 4,
                         spec.liquid_material_fingerprint, {injections, 1},
                         tcr.snapshot(), 1U << 20))
    return false;
  const auto save = [](RestartCellRecordsView records, std::uint64_t step) {
    RestartImage image;
    image.source_format_version = 5;
    image.step = step;
    image.backward_euler_recovery = false;
    image.cell_record_identity = records.identity;
    image.cell_record_bytes = records.record_bytes;
    if (records.values.size)
      image.cell_records.assign(records.values.data,
                                records.values.data + records.values.size);
    image.cell_record_lengths.assign(records.variable_cell_bytes.data,
                                     records.variable_cell_bytes.data +
                                         records.variable_cell_bytes.size);
    return image;
  };
  const auto initial = save(history.snapshot(), 0);
  if (!injector.begin_trial(0, 1).succeeded())
    return false;
  std::array<hundun::v04::detail::ProductSprayHistory::Parcel, 2> parcels;
  if (!injector.candidate_at(0, parcels[0].parcel) ||
      !injector.candidate_at(1, parcels[1].parcel))
    return false;
  parcels[0].parcel.owner_global_cell = 0;
  parcels[0].tab_deformation = -.25;
  parcels[0].tab_deformation_rate_per_s = 12;
  parcels[0].breakup_ordinal = UINT64_C(9007199254740997);
  parcels[1].breakup_ordinal = UINT64_MAX - 1;
  if (!stage_tcr())
    return false;
  hot_allocation_count = 0;
  count_hot_allocations = true;
  const auto staged = history.stage_next({parcels.data(), parcels.size()}, 0,
                                         tcr.prepared_snapshot());
  const auto ready = history.preflight_commit();
  count_hot_allocations = false;
  auto accepted = save(history.snapshot(), 0);
  bool passed = expect(staged && ready && hot_allocation_count == 0 &&
                           history.accepted_parcels().size == 0 &&
                           history.prepared_parcels().size == 2 &&
                           accepted.cell_records == initial.cell_records &&
                           injector.committed_state() == original,
                       "native parcel, encoded history and injector counters "
                       "remain pending together");
  if (!passed)
    return false;
  history.discard();
  tcr.discard();
  passed &= expect(!history.preflight_commit() && !injector.trial_active() &&
                       injector.committed_state() == original,
                   "gas rejection withdraws all pending spray participants");
  if (!injector.begin_trial(0, 1).succeeded() || !stage_tcr() ||
      !history.stage_next({parcels.data(), parcels.size()}, 0,
                          tcr.prepared_snapshot()) ||
      !history.preflight_commit())
    return false;
  history.commit();
  tcr.commit();
  auto image = save(history.snapshot(), 1);
  const auto emitted = injector.committed_state();
  passed &= expect(history.accepted_parcels().size == 2 &&
                       emitted.next_ordinal == original.next_ordinal + 2 &&
                       image.cell_record_lengths[0] ==
                           24 + 144 + (coupled_tcr ? 120 : 0) &&
                       image.cell_record_lengths[7] ==
                           24 + 144 + 24 + (coupled_tcr ? 120 : 0),
                   "one publication advances parcels and exact counters into "
                   "cell-partitioned V5 records");
  if (!expect(bool(history.stage_restore(initial)),
              "initial complete-history snapshot stages"))
    return false;
  if (!tcr.stage_restore_records(history.prepared_tcr_records(), 0))
    return false;
  history.commit();
  tcr.commit();
  auto corrupt = image;
  corrupt.cell_records[0] = 2; // embedded cell step disagrees with native clock
  passed &= expect(
      !history.stage_restore(corrupt) && history.accepted_parcels().size == 0 &&
          injector.committed_state() == original,
      "corrupt cell clock rejects before any parcel or injector is restored");
  hot_allocation_count = 0;
  count_hot_allocations = true;
  const auto restore = history.stage_restore(image);
  const auto restore_ready = history.preflight_commit();
  const auto tcr_restore =
      tcr.stage_restore_records(history.prepared_tcr_records(), 1);
  count_hot_allocations = false;
  passed &= expect(restore && restore_ready && tcr_restore &&
                       hot_allocation_count == 0 &&
                       injector.committed_state() == original &&
                       history.accepted_parcels().size == 0,
                   "valid V5 restore stages all spray state without allocation "
                   "or early publication");
  if (!restore || !restore_ready || !tcr_restore)
    return false;
  history.commit();
  tcr.commit();
  const auto restored = history.accepted_parcels();
  passed &= expect(
      restored.size == 2 &&
          same_parcel(restored.data[0].parcel, parcels[0].parcel) &&
          same_parcel(restored.data[1].parcel, parcels[1].parcel) &&
          restored.data[0].tab_deformation == -.25 &&
          restored.data[0].tab_deformation_rate_per_s == 12 &&
          restored.data[0].breakup_ordinal == UINT64_C(9007199254740997) &&
          restored.data[1].breakup_ordinal == UINT64_MAX - 1 &&
          injector.committed_state() == emitted,
      "V5 restore preserves full parcel IDs, TAB and uint64 breakup/injection "
      "lineage");
  if (coupled_tcr) {
    auto damaged = image;
    damaged.cell_records[24 + 64] ^=
        1; // valid outer record, wrong TCR mapping identity
    const auto parcel_restore = history.stage_restore(damaged);
    const auto tcr_restore =
        tcr.stage_restore_records(history.prepared_tcr_records(), 1);
    history.discard();
    tcr.discard();
    const auto unchanged = history.snapshot();
    passed &=
        expect(parcel_restore && !tcr_restore && !injector.trial_active() &&
                   injector.committed_state() == emitted &&
                   unchanged.values.size == image.cell_records.size() &&
                   std::equal(image.cell_records.begin(),
                              image.cell_records.end(), unchanged.values.data),
               "late typed TCR rejection withdraws prepared injector/parcel "
               "restore without changing accepted bytes");
  }
  if (coupled_tcr)
    passed &= expect(tcr.accepted(7).revision.accepted_step == 1 &&
                         tcr.accepted(7).branch_sign == -1 &&
                         std::abs(tcr.accepted(7).control - 1.0 / 3) < 1e-14,
                     "combined V5 payload restores typed TCR branch history "
                     "with spray state");
  return passed;
}

bool test_staged_injector_restore() {
  DeterministicInjector injector, reference;
  const auto spec = base_injector();
  const InjectorCommittedState original{.25, 100};
  const InjectorCommittedState restored{.375, UINT64_C(9007199254740993)};
  if (!injector.reserve(4) || !injector.configure(spec, original) ||
      !reference.reserve(4) || !reference.configure(spec, restored))
    return false;
  InjectorCommittedState candidate;
  hot_allocation_count = 0;
  count_hot_allocations = true;
  const auto staged = injector.stage_restore(restored);
  const auto prepared = injector.prepared_state(candidate);
  const auto preflight = injector.preflight_commit();
  count_hot_allocations = false;
  bool passed =
      expect(staged && prepared && preflight && candidate == restored &&
                 injector.committed_state() == original &&
                 injector.candidate_count() == 0 && hot_allocation_count == 0,
             "restore counters stage without injecting parcels or publishing "
             "accepted state");
  if (!passed)
    return false;
  passed &= expect(
      injector.rollback_trial() && injector.committed_state() == original &&
          !injector.prepared_state(candidate),
      "failed native restore can withdraw prepared injector counters");
  passed &=
      expect(!injector.stage_restore({-1, 7}) && !injector.trial_active() &&
                 injector.committed_state() == original,
             "invalid restart remainder leaves injector untouched");
  passed &=
      expect(injector.stage_restore(restored) && injector.commit_trial() &&
                 injector.committed_state() == restored,
             "restart publication preserves counters above 2^53 exactly");
  const auto a = injector.begin_trial(42, 1), b = reference.begin_trial(42, 1);
  passed &= expect(a.succeeded() && b.succeeded() && same_report(a, b),
                   "restored counters produce the same next injection report");
  for (std::size_t i = 0; i < a.parcel_count; ++i) {
    SprayParcelState x, y;
    passed &=
        expect(injector.candidate_at(i, x) && reference.candidate_at(i, y) &&
                   same_parcel(x, y),
               "post-restore injector emits exact reference parcel identities");
  }
  return passed;
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
  ParcelContainer container;
  DeterministicInjector injector;
  passed &= expect(static_cast<bool>(container.reserve(2)) &&
      static_cast<bool>(container.begin_trial()) &&
      static_cast<bool>(container.stage_add(parcel(1,2,0))) &&
      static_cast<bool>(container.stage_tab_state({1,2},0.25,12.0)) &&
      static_cast<bool>(container.commit_trial()) &&
      static_cast<bool>(injector.reserve(3)) &&
      static_cast<bool>(injector.configure(base_injector(),{0.25,9})),
      "accepted lifecycle fixture");
  const DeterministicInjector* injectors[]{&injector};
  const std::uint64_t ordinals[]{7};
  const hundun::v04::portable::Revision revision{12,31,1};
  const auto snapshot=export_parcel_lifecycle(container,{injectors,1},{ordinals,1},
      revision,0.125,43);
  const auto restore=prepare_parcel_lifecycle_restore(snapshot.values,revision,2,1);
  passed &= expect(snapshot.available && restore.available &&
      restore.values.parcels.tab_deformation[0]==0.25 &&
      restore.values.parcels.tab_deformation_rate_per_s[0]==12.0 &&
      restore.values.breakup_ordinals[0]==7 &&
      restore.values.injectors[0].accepted==InjectorCommittedState{0.25,9} &&
      restore.values.parcel_rng_seed==43 && restore.values.accepted_time_s==0.125,
      "owning restore candidate covers TAB, breakup, injector and accepted RNG clocks");
  // Live trial values are never exported. Restoring accepted injector counters
  // reproduces the same next accepted-step candidate IDs and residual exactly.
  passed &= expect(static_cast<bool>(container.begin_trial()) &&
      static_cast<bool>(container.stage_tab_state({1,2},9.0,99.0)), "unaccepted TAB trial");
  const auto original_injection=injector.begin_trial(13,1.0);
  const auto during_trial=export_parcel_lifecycle(container,{injectors,1},{ordinals,1},
      revision,0.125,43);
  passed &= expect(original_injection.succeeded() && during_trial.available &&
      same_snapshot_bytes(during_trial.values.parcels,snapshot.values.parcels) &&
      during_trial.values.injectors[0].accepted==snapshot.values.injectors[0].accepted,
      "active trial TAB and injection candidates do not enter accepted snapshots");
  DeterministicInjector restored_injector;
  passed &= expect(static_cast<bool>(restored_injector.reserve(3)) &&
      static_cast<bool>(restored_injector.configure(restore.values.injectors[0].spec,
          restore.values.injectors[0].accepted)) &&
      restored_injector.begin_trial(13,1.0).succeeded(), "restored injector prepares retry");
  for(std::size_t i=0;i<injector.candidate_count();++i) {
    SprayParcelState a,b;
    passed &= expect(injector.candidate_at(i,a) && restored_injector.candidate_at(i,b) &&
        same_parcel(a,b), "continuous and restored injection candidate identities match exactly");
  }
  passed &= expect(static_cast<bool>(container.rollback_trial()) &&
      static_cast<bool>(injector.rollback_trial()), "unaccepted trial rollback");
  auto corrupt=snapshot.values;
  corrupt.parcels.tab_deformation[0]=std::numeric_limits<double>::quiet_NaN();
  passed &= expect(!prepare_parcel_lifecycle_restore(corrupt,revision,2,1).available,
      "nonfinite TAB history rejects the entire restore candidate");
  corrupt=snapshot.values;
  corrupt.injectors[0].accepted.residual_mass_kg=-1;
  passed &= expect(!prepare_parcel_lifecycle_restore(corrupt,revision,2,1).available,
      "invalid injector remainder rejects the entire restore candidate");
  corrupt=snapshot.values;
  corrupt.parcel_rng_algorithm_version=2;
  passed &= expect(!prepare_parcel_lifecycle_restore(corrupt,revision,2,1).available &&
      !prepare_parcel_lifecycle_restore(snapshot.values,{12,32,1},2,1).available &&
      !prepare_parcel_lifecycle_restore(snapshot.values,revision,0,1).available &&
      !prepare_parcel_lifecycle_restore(snapshot.values,revision,2,0).available,
      "unknown RNG, stale revision and capacity shortages cannot publish partial restore");
  passed &= test_injector_residual_retry_and_identity();
  passed &= test_staged_injector_restore();
  passed &= test_native_spray_history(false);
  passed &= test_native_spray_history(true);
  passed &= test_micro_mass_residual_scale();
  passed &= test_cone_zero_flow_and_capacity_failure();
  return passed ? 0 : 1;
}

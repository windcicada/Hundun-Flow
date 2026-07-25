// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "Checkpoint v2 test access is unavailable when tests are disabled"
#endif

#include "checkpoint_v2_detail.hpp"
#include "../../runtime/src/field_epoch_test_access.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <vector>

namespace hundun::runtime::checkpoint_v2 {
struct Manifest;
}

namespace hundun::flow::test {

enum class CheckpointV2PreparationPoint : std::uint8_t {
  none,
  local_layout,
  local_topology,
  local_geometry,
  topology_common,
  geometry_common,
  resolved_case,
  boundary_registry,
  field_schema,
  common_authority,
  final_success_boundary,
  path,
  rank_path
};

void set_checkpoint_v2_preparation_fault(
    CheckpointV2PreparationPoint, std::uint32_t calls_before = 0U) noexcept;
void set_ideal_gas_restore_snapshot_preparation_fault(int rank) noexcept;

struct CheckpointV2PathCodeObservation final {
  std::uint64_t success_code{};
  std::uint64_t candidate_code{};
  bool decoded{};
  bool success{};
  int rank{-1};
  int reason{-1};
};

CheckpointV2PathCodeObservation
checkpoint_v2_path_code_observation_for_test(int size, int rank, int reason,
                                             bool local_success,
                                             std::uint64_t selected) noexcept;
std::vector<std::uint8_t> checkpoint_v2_encode_global_payload_for_test(
    AcceptedStepMetadata, const TimeControlState &,
    const std::optional<IdealGasClosureState> &);
std::vector<std::uint8_t>
checkpoint_v2_encode_rank_payload_for_test(const FlowState &,
                                           std::uint64_t &logical_bytes);
bool checkpoint_v2_authenticate_global_payload_for_test(
    const std::vector<std::uint8_t> &, std::uint64_t expected_crc,
    std::uint64_t expected_size) noexcept;
bool checkpoint_v2_authenticate_rank_payload_for_test(
    const std::vector<std::uint8_t> &, const FlowState &,
    std::uint64_t expected_crc) noexcept;
bool checkpoint_v2_authenticate_rank_wrapper_for_test(
    const std::vector<std::uint8_t> &, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::int32_t expected_rank,
    std::int32_t expected_rank_count,
    std::uint64_t expected_payload_size) noexcept;
bool checkpoint_v2_authenticate_manifest_for_test(
    const std::vector<std::uint8_t> &, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size,
    const runtime::checkpoint_v2::Manifest &expected) noexcept;
bool checkpoint_v2_authenticate_manifest_limits_for_test(
    const std::vector<std::uint8_t> &, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::uint64_t expected_rank_count,
    std::uint64_t expected_global_payload_size) noexcept;
bool checkpoint_v2_authenticate_completed_marker_for_test(
    const std::vector<std::uint8_t> &,
    std::uint64_t expected_manifest_actual_size,
    std::uint64_t expected_manifest_crc64,
    std::uint64_t expected_common_fingerprint) noexcept;
std::uint64_t
checkpoint_v2_field_schema_fingerprint_for_test(const runtime::FieldRegistry &,
                                                const FlowFieldIds &);

struct CheckpointV2DeepSnapshot final {
  std::array<std::vector<std::vector<std::uint64_t>>, 4> storage_bits;
  std::array<std::uint64_t, 4> generations{};
  AcceptedStepMetadata metadata{};
  bool rollback_snapshot_valid{};
  bool attempt_active{};
  std::uint64_t attempt_identity{};
  std::uint64_t diagnostic_mutation_identity{};
  bool commit_prepared{};
  AcceptedStepMetadata prepared_metadata{};
};

inline bool checkpoint_v2_deep_snapshot_equal(
    const CheckpointV2DeepSnapshot &left,
    const CheckpointV2DeepSnapshot &right) noexcept {
  const auto metadata_equal = [](AcceptedStepMetadata lhs,
                                 AcceptedStepMetadata rhs) noexcept {
    const auto bits = [](double value) noexcept {
      std::uint64_t result{};
      std::memcpy(&result, &value, sizeof(result));
      return result;
    };
    return lhs.step == rhs.step && lhs.order == rhs.order &&
           bits(lhs.time_s) == bits(rhs.time_s) &&
           bits(lhs.dt_s) == bits(rhs.dt_s) &&
           bits(lhs.previous_dt_s) == bits(rhs.previous_dt_s);
  };
  return left.storage_bits == right.storage_bits &&
         left.generations == right.generations &&
         metadata_equal(left.metadata, right.metadata) &&
         left.rollback_snapshot_valid == right.rollback_snapshot_valid &&
         left.attempt_active == right.attempt_active &&
         left.attempt_identity == right.attempt_identity &&
         left.diagnostic_mutation_identity ==
             right.diagnostic_mutation_identity &&
         left.commit_prepared == right.commit_prepared &&
         metadata_equal(left.prepared_metadata, right.prepared_metadata);
}

inline bool checkpoint_v2_failed_read_preserved_values(
    const CheckpointV2DeepSnapshot &before,
    const CheckpointV2DeepSnapshot &after) noexcept {
  auto normalized = after;
  normalized.generations = before.generations;
  normalized.diagnostic_mutation_identity =
      before.diagnostic_mutation_identity;
  if (!checkpoint_v2_deep_snapshot_equal(before, normalized))
    return false;
  for (std::size_t index = 0; index < before.generations.size(); ++index) {
    if (before.generations[index] ==
            std::numeric_limits<std::uint64_t>::max() ||
        after.generations[index] != before.generations[index] + 1U)
      return false;
  }
  return before.diagnostic_mutation_identity !=
             std::numeric_limits<std::uint64_t>::max() &&
         after.diagnostic_mutation_identity ==
             before.diagnostic_mutation_identity + 1U;
}

class CheckpointV2TestAccess final {
public:
  static CheckpointV2DeepSnapshot snapshot(const FlowState &state) {
    if (!state.impl_)
      throw runtime::Error("Checkpoint v2 test FlowState has been moved from");
    CheckpointV2DeepSnapshot result;
    const std::array<const runtime::FieldStorage *, 4> storages{
        &state.impl_->history, &state.impl_->committed, &state.impl_->trial,
        &state.impl_->rollback_snapshot};
    const std::vector<runtime::FieldId> fields{
        state.impl_->fields.density,
        state.impl_->fields.velocity,
        state.impl_->fields.mechanical_pressure,
        state.impl_->fields.face_velocity,
        state.impl_->fields.face_mass_flux};
    std::vector<runtime::FieldId> ordered(fields);
    ordered.insert(ordered.end(),
                   state.impl_->fields.transported_cell_fields.begin(),
                   state.impl_->fields.transported_cell_fields.end());
    runtime::FieldAccessPlan read_plan(*state.impl_->registry);
    constexpr runtime::PhaseId phase = 2300U;
    constexpr runtime::ActorId actor = 2300U;
    for (const auto field : ordered)
      read_plan.declare_access(phase, actor, field,
                               runtime::AccessMode::read);
    read_plan.freeze();

    for (std::size_t layer = 0; layer < storages.size(); ++layer) {
      const auto &storage = *storages[layer];
      result.generations[layer] =
          runtime::detail::FieldEpochTestAccess::generation(storage);
      for (const auto field : ordered) {
        const auto &descriptor = state.impl_->registry->descriptor(field);
        std::vector<std::uint64_t> values;
        if (descriptor.space == runtime::FunctionSpace::cell_average) {
          const auto view = storage.view<double>(field);
          const auto extent = view.interior_extent();
          const int ghost = view.ghost_width();
          const auto nx = static_cast<std::size_t>(extent.x + 2 * ghost);
          const auto ny = static_cast<std::size_t>(extent.y + 2 * ghost);
          const auto nz = static_cast<std::size_t>(extent.z + 2 * ghost);
          values.reserve(nx * ny * nz * descriptor.components);
          for (int k = -ghost; k < extent.z + ghost; ++k)
            for (int j = -ghost; j < extent.y + ghost; ++j)
              for (int i = -ghost; i < extent.x + ghost; ++i)
                for (std::uint32_t component = 0;
                     component < descriptor.components; ++component)
                  values.push_back(bits(view(
                      i, j, k, static_cast<int>(component))));
        } else if (descriptor.space == runtime::FunctionSpace::face_value) {
          const auto view = storage.acquire_face_read<double>(
              read_plan, phase, actor, field);
          values.reserve(view.face_count() * descriptor.components);
          for (std::size_t face = 0; face < view.face_count(); ++face)
            for (std::uint32_t component = 0;
                 component < descriptor.components; ++component)
              values.push_back(
                  bits(view(face, static_cast<int>(component))));
        } else {
          throw runtime::Error(
              "Checkpoint v2 test encountered unsupported field space");
        }
        result.storage_bits[layer].push_back(std::move(values));
      }
    }
    result.metadata = state.impl_->metadata;
    result.rollback_snapshot_valid = state.impl_->rollback_snapshot_valid;
    result.attempt_active = state.impl_->attempt_active;
    result.attempt_identity = state.impl_->attempt_identity;
    result.diagnostic_mutation_identity =
        state.impl_->diagnostic_mutation_identity;
    result.commit_prepared = state.impl_->commit_prepared;
    result.prepared_metadata = state.impl_->prepared_metadata;
    return result;
  }

  static void set_committed_density_ghost(FlowState &state, double value) {
    auto view =
        state.impl_->committed.view<double>(state.impl_->fields.density);
    view(-1, 0, 0, 0) = value;
  }

  static void set_rollback_density_ghost(FlowState &state, double value) {
    auto view = state.impl_->rollback_snapshot.view<double>(
        state.impl_->fields.density);
    view(-1, 0, 0, 0) = value;
  }

  static void force_generation(FlowState &state, std::size_t layer,
                               std::uint64_t generation) {
    auto storages = mutable_storages(state);
    if (layer >= storages.size())
      throw runtime::Error("Checkpoint v2 test layer is out of range");
    runtime::detail::FieldEpochTestAccess::force_generation(*storages[layer],
                                                            generation);
  }

  static std::array<runtime::FieldView<const double>, 4>
  density_views(const FlowState &state) {
    if (!state.impl_)
      throw runtime::Error("Checkpoint v2 test FlowState has been moved from");
    const auto field = state.impl_->fields.density;
    const runtime::FieldStorage &history = state.impl_->history;
    const runtime::FieldStorage &committed = state.impl_->committed;
    const runtime::FieldStorage &trial = state.impl_->trial;
    const runtime::FieldStorage &rollback = state.impl_->rollback_snapshot;
    return {history.view<double>(field), committed.view<double>(field),
            trial.view<double>(field), rollback.view<double>(field)};
  }

  static bool all_cell_ghosts_are_positive_zero(const FlowState &state) {
    if (!state.impl_)
      throw runtime::Error("Checkpoint v2 test FlowState has been moved from");
    const std::array<const runtime::FieldStorage *, 4> storages{
        &state.impl_->history, &state.impl_->committed, &state.impl_->trial,
        &state.impl_->rollback_snapshot};
    std::vector<runtime::FieldId> fields{
        state.impl_->fields.density, state.impl_->fields.velocity,
        state.impl_->fields.mechanical_pressure};
    fields.insert(fields.end(),
                  state.impl_->fields.transported_cell_fields.begin(),
                  state.impl_->fields.transported_cell_fields.end());
    for (const auto *storage : storages)
      for (const auto field : fields) {
        const auto view = storage->view<double>(field);
        const auto extent = view.interior_extent();
        const int ghost = view.ghost_width();
        const auto components =
            state.impl_->registry->descriptor(field).components;
        for (int k = -ghost; k < extent.z + ghost; ++k)
          for (int j = -ghost; j < extent.y + ghost; ++j)
            for (int i = -ghost; i < extent.x + ghost; ++i) {
              const bool interior = i >= 0 && i < extent.x && j >= 0 &&
                                    j < extent.y && k >= 0 && k < extent.z;
              if (interior)
                continue;
              for (std::uint32_t component = 0; component < components;
                   ++component)
                if (bits(view(i, j, k, static_cast<int>(component))) != 0U)
                  return false;
            }
      }
    return true;
  }

private:
  static std::uint64_t bits(double value) noexcept {
    std::uint64_t result{};
    std::memcpy(&result, &value, sizeof(result));
    return result;
  }

  static std::array<runtime::FieldStorage *, 4>
  mutable_storages(FlowState &state) {
    if (!state.impl_)
      throw runtime::Error("Checkpoint v2 test FlowState has been moved from");
    return {&state.impl_->history, &state.impl_->committed,
            &state.impl_->trial, &state.impl_->rollback_snapshot};
  }
};

inline bool checkpoint_v2_deep_snapshot_oracle_is_mutation_sensitive(
    FlowState &state) {
  const auto baseline = CheckpointV2TestAccess::snapshot(state);
  auto ordinary = baseline;
  ordinary.storage_bits[1][0][0] ^= 1U;
  auto nested = baseline;
  nested.storage_bits[1].back().back() ^= 1U;
  auto ghost = baseline;
  ghost.storage_bits[1][0].front() ^= 1U;
  auto rollback = baseline;
  rollback.storage_bits[3][0].front() ^= 1U;
  return checkpoint_v2_deep_snapshot_equal(baseline, baseline) &&
         !checkpoint_v2_deep_snapshot_equal(baseline, ordinary) &&
         !checkpoint_v2_deep_snapshot_equal(baseline, nested) &&
         !checkpoint_v2_deep_snapshot_equal(baseline, ghost) &&
         !checkpoint_v2_deep_snapshot_equal(baseline, rollback);
}

}

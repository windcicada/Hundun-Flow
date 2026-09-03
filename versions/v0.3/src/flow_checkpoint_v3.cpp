// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_checkpoint_v3.hpp"

#include "flow_checkpoint_v3_detail.hpp"
#include "flow_checkpoint_v2_detail.hpp"
#include "flow_density_closure_detail.hpp"
#include "rt_checkpoint_v2_protocol_detail.hpp"

#include "hundun/cfg_resolved_case_v3_loader.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr std::array<std::uint8_t, 16> kManifestMagic{
    'H', 'U', 'N', 'D', 'U', 'N', '-', 'C',
    'K', 'P', 'T', '-', 'V', '3', 0,   0};
constexpr std::uint32_t kManifestSchema = 3U;
constexpr std::uint64_t kMaximumAuthorityRows = UINT64_C(100000000);
constexpr std::uint32_t kWaleTransientSchemaVersion = 1U;
constexpr std::uint64_t kWaleSectionBytes = 47U;
constexpr std::uint64_t kIdealGasSectionBytes = 30U;

std::atomic<std::uint64_t> checkpoint_logical_io_bytes{};

void record_checkpoint_logical_bytes(std::uint64_t bytes) {
  auto value = checkpoint_logical_io_bytes.load(std::memory_order_relaxed);
  for (;;) {
    if (bytes > std::numeric_limits<std::uint64_t>::max() - value)
      throw runtime::Error("Checkpoint v3 logical I/O counter would overflow");
    if (checkpoint_logical_io_bytes.compare_exchange_weak(
            value, value + bytes, std::memory_order_relaxed,
            std::memory_order_relaxed))
      return;
  }
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(runtime::Box3 left, runtime::Box3 right) noexcept {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

bool valid_rows(
    bool available,
    const std::vector<detail::CheckpointV3AuthorityGradient> &rows) noexcept {
  if (!available) {
    return rows.empty();
  }
  if (rows.size() > static_cast<std::size_t>(kMaximumAuthorityRows)) {
    return false;
  }
  for (std::size_t index = 0U; index < rows.size(); ++index) {
    if (!std::isfinite(rows[index].value) ||
        (index != 0U && rows[index - 1U].link >= rows[index].link)) {
      return false;
    }
  }
  return true;
}

bool same_layout(
    const std::vector<detail::CheckpointV3AuthorityGradient> &left,
    const std::vector<detail::CheckpointV3AuthorityGradient> &right) noexcept;

bool valid_rank_authority(
    const detail::CheckpointV3RankAuthority &rank) noexcept {
  return rank.rank >= 0 && rank.owned_box.begin.x >= 0 &&
         rank.owned_box.begin.y >= 0 && rank.owned_box.begin.z >= 0 &&
         rank.owned_box.end.x > rank.owned_box.begin.x &&
         rank.owned_box.end.y > rank.owned_box.begin.y &&
         rank.owned_box.end.z > rank.owned_box.begin.z &&
         !rank.state_filename.empty() && rank.state_actual_bytes != 0U &&
         valid_rows(rank.history_available, rank.history) &&
         valid_rows(rank.committed_available, rank.committed) &&
         (!rank.history_available || !rank.committed_available ||
          same_layout(rank.history, rank.committed));
}

bool same_layout(
    const std::vector<detail::CheckpointV3AuthorityGradient> &left,
    const std::vector<detail::CheckpointV3AuthorityGradient> &right) noexcept {
  if (left.size() != right.size()) {
    return false;
  }
  for (std::size_t index = 0U; index < left.size(); ++index) {
    if (left[index].link != right[index].link) {
      return false;
    }
  }
  return true;
}

bool process_grid_matches(runtime::Int3 grid, std::int32_t ranks) noexcept {
  if (grid.x <= 0 || grid.y <= 0 || grid.z <= 0 || ranks <= 0) {
    return false;
  }
  const auto x = static_cast<std::uint64_t>(grid.x);
  const auto y = static_cast<std::uint64_t>(grid.y);
  const auto z = static_cast<std::uint64_t>(grid.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y) {
    return false;
  }
  const auto xy = x * y;
  if (xy > std::numeric_limits<std::uint64_t>::max() / z) {
    return false;
  }
  return xy * z == static_cast<std::uint64_t>(ranks);
}

bool known_presence(CheckpointV3Presence presence) noexcept {
  switch (presence) {
  case CheckpointV3Presence::constant_static_ibm:
  case CheckpointV3Presence::constant_body_fitted_wale:
  case CheckpointV3Presence::constant_static_ibm_wale:
  case CheckpointV3Presence::material_static_ibm:
  case CheckpointV3Presence::material_body_fitted_wale:
  case CheckpointV3Presence::material_static_ibm_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm:
  case CheckpointV3Presence::ideal_gas_body_fitted_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm_wale:
    return true;
  }
  return false;
}

bool profile_has_ibm(CheckpointV3Presence presence) noexcept {
  switch (presence) {
  case CheckpointV3Presence::constant_static_ibm:
  case CheckpointV3Presence::constant_static_ibm_wale:
  case CheckpointV3Presence::material_static_ibm:
  case CheckpointV3Presence::material_static_ibm_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm:
  case CheckpointV3Presence::ideal_gas_static_ibm_wale:
    return true;
  case CheckpointV3Presence::constant_body_fitted_wale:
  case CheckpointV3Presence::material_body_fitted_wale:
  case CheckpointV3Presence::ideal_gas_body_fitted_wale:
    return false;
  }
  return false;
}

bool profile_has_wale(CheckpointV3Presence presence) noexcept {
  switch (presence) {
  case CheckpointV3Presence::constant_body_fitted_wale:
  case CheckpointV3Presence::constant_static_ibm_wale:
  case CheckpointV3Presence::material_body_fitted_wale:
  case CheckpointV3Presence::material_static_ibm_wale:
  case CheckpointV3Presence::ideal_gas_body_fitted_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm_wale:
    return true;
  case CheckpointV3Presence::constant_static_ibm:
  case CheckpointV3Presence::material_static_ibm:
  case CheckpointV3Presence::ideal_gas_static_ibm:
    return false;
  }
  return false;
}

bool profile_has_material_fields(CheckpointV3Presence presence) noexcept {
  return static_cast<std::uint8_t>(presence) >=
             static_cast<std::uint8_t>(
                 CheckpointV3Presence::material_static_ibm) &&
         known_presence(presence);
}

bool profile_has_ideal_gas(CheckpointV3Presence presence) noexcept {
  return static_cast<std::uint8_t>(presence) >=
             static_cast<std::uint8_t>(
                 CheckpointV3Presence::ideal_gas_static_ibm) &&
         known_presence(presence);
}

bool valid_ideal_gas_state(const IdealGasClosureState &state) noexcept {
  if (!(state.thermodynamic_pressure_pa > 0.0) ||
      !std::isfinite(state.thermodynamic_pressure_pa) ||
      state.revision == std::numeric_limits<std::uint64_t>::max()) {
    return false;
  }
  switch (state.mode) {
  case IdealGasPressureMode::open_fixed:
    return !state.target_mass_kg.has_value();
  case IdealGasPressureMode::closed_dynamic:
    return state.target_mass_kg.has_value() &&
           *state.target_mass_kg > 0.0 &&
           std::isfinite(*state.target_mass_kg);
  }
  return false;
}

bool canonical_absent_ibm(
    const detail::CheckpointV3RankAuthority &rank) noexcept {
  return std::all_of(rank.fingerprints.begin(), rank.fingerprints.end(),
                     [](std::uint64_t value) { return value == 0U; }) &&
         !rank.history_available && !rank.committed_available &&
         rank.history.empty() && rank.committed.empty();
}

bool valid_wale_identity(
    const detail::CheckpointV3WaleIdentity &identity) noexcept {
  return std::isfinite(identity.coefficient) &&
         identity.coefficient >= 1.0e-6 && identity.coefficient <= 1.0 &&
         std::isfinite(identity.turbulent_prandtl) &&
         identity.turbulent_prandtl >= 0.1 &&
         identity.turbulent_prandtl <= 10.0 &&
         std::isfinite(identity.turbulent_schmidt) &&
         identity.turbulent_schmidt >= 0.1 &&
         identity.turbulent_schmidt <= 10.0 &&
         identity.transient_schema_version ==
             kWaleTransientSchemaVersion &&
         identity.transient_fields ==
             std::array<detail::CheckpointV3WaleTransientField, 3>{
                 detail::CheckpointV3WaleTransientField::nu_t_m2_per_s,
                 detail::CheckpointV3WaleTransientField::mu_sgs_pa_s,
                 detail::CheckpointV3WaleTransientField::mu_eff_pa_s};
}

bool add_section_bytes(std::uint64_t &total, std::uint64_t bytes) noexcept {
  if (bytes > std::numeric_limits<std::uint64_t>::max() - total) {
    return false;
  }
  total += bytes;
  return true;
}

std::optional<std::uint64_t> ibm_section_bytes(
    const std::vector<detail::CheckpointV3RankAuthority> &ranks) noexcept {
  std::uint64_t total = 8U;
  for (const auto &rank : ranks) {
    constexpr std::uint64_t fixed_rank_bytes = 94U;
    if (rank.history.size() >
            std::numeric_limits<std::uint64_t>::max() / 16U ||
        rank.committed.size() >
            std::numeric_limits<std::uint64_t>::max() / 16U ||
        !add_section_bytes(total, fixed_rank_bytes) ||
        !add_section_bytes(
            total, 16U * static_cast<std::uint64_t>(rank.history.size())) ||
        !add_section_bytes(
            total,
            16U * static_cast<std::uint64_t>(rank.committed.size()))) {
      return std::nullopt;
    }
  }
  return total;
}

bool valid_manifest(const detail::CheckpointV3Manifest &value) noexcept {
  if (!known_presence(value.presence) ||
      !process_grid_matches(value.process_grid, value.rank_count) ||
      !std::isfinite(value.metadata.time_s) || value.metadata.time_s < 0.0 ||
      !std::isfinite(value.metadata.dt_s) || value.metadata.dt_s <= 0.0 ||
      !std::isfinite(value.metadata.previous_dt_s) ||
      value.metadata.previous_dt_s < 0.0 ||
      !std::isfinite(value.control.proposed_next_dt_s) ||
      !(value.control.proposed_next_dt_s > 0.0) ||
      value.control.last_retry_count > 8U ||
      value.ranks.size() != static_cast<std::size_t>(value.rank_count)) {
    return false;
  }
  for (std::size_t index = 0U; index < value.ranks.size(); ++index) {
    const auto &rank = value.ranks[index];
    if (rank.rank != static_cast<std::int32_t>(index) ||
        !valid_rank_authority(rank)) {
      return false;
    }
  }
  switch (value.presence) {
  case CheckpointV3Presence::constant_static_ibm:
    if (value.ibm_section_count != 0U || value.ibm_section_bytes != 0U ||
        value.wale_section_count != 0U || value.wale_section_bytes != 0U ||
        value.wale.has_value() || value.ideal_gas_section_count != 0U ||
        value.ideal_gas_section_bytes != 0U || value.ideal_gas.has_value()) {
      return false;
    }
    break;
  case CheckpointV3Presence::constant_body_fitted_wale:
  case CheckpointV3Presence::constant_static_ibm_wale:
  case CheckpointV3Presence::material_static_ibm:
  case CheckpointV3Presence::material_body_fitted_wale:
  case CheckpointV3Presence::material_static_ibm_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm:
  case CheckpointV3Presence::ideal_gas_body_fitted_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm_wale:
    if (profile_has_ibm(value.presence)) {
      const auto expected_ibm_bytes = ibm_section_bytes(value.ranks);
      if (value.ibm_section_count != 1U || !expected_ibm_bytes.has_value() ||
          value.ibm_section_bytes != *expected_ibm_bytes) {
        return false;
      }
    } else if (value.ibm_section_count != 0U ||
               value.ibm_section_bytes != 0U ||
               !std::all_of(value.ranks.begin(), value.ranks.end(),
                            canonical_absent_ibm)) {
      return false;
    }
    if (profile_has_wale(value.presence)) {
      if (value.wale_section_count != 1U ||
          value.wale_section_bytes != kWaleSectionBytes ||
          !value.wale.has_value() || !valid_wale_identity(*value.wale)) {
        return false;
      }
    } else if (value.wale_section_count != 0U ||
               value.wale_section_bytes != 0U || value.wale.has_value()) {
      return false;
    }
    if (profile_has_ideal_gas(value.presence)) {
      if (value.ideal_gas_section_count != 1U ||
          value.ideal_gas_section_bytes != kIdealGasSectionBytes ||
          !value.ideal_gas.has_value() ||
          !valid_ideal_gas_state(*value.ideal_gas)) {
        return false;
      }
    } else if (value.ideal_gas_section_count != 0U ||
               value.ideal_gas_section_bytes != 0U ||
               value.ideal_gas.has_value()) {
      return false;
    }
    break;
  }
  if (value.metadata.step == 0U) {
    if (value.metadata.time_s != 0.0 ||
        value.metadata.previous_dt_s != 0.0 ||
        value.metadata.order != MomentumTimeOrder::backward_euler) {
      return false;
    }
  } else {
    if (!(value.metadata.time_s > 0.0) ||
        !(value.metadata.previous_dt_s > 0.0)) {
      return false;
    }
    if ((value.metadata.step == 1U &&
         value.metadata.order != MomentumTimeOrder::backward_euler) ||
        (value.metadata.step > 1U &&
         value.metadata.order != MomentumTimeOrder::bdf2)) {
      return false;
    }
  }
  return true;
}

void require_valid(const detail::CheckpointV3Manifest &value) {
  if (!valid_manifest(value)) {
    throw runtime::Error("Checkpoint v3 manifest is invalid");
  }
}

void append_rows(
    runtime::checkpoint_v2::Encoder &encoder,
    const std::vector<detail::CheckpointV3AuthorityGradient> &rows) {
  encoder.u64(static_cast<std::uint64_t>(rows.size()));
  for (const auto &row : rows) {
    encoder.u64(row.link);
    encoder.f64(row.value);
  }
}

std::vector<detail::CheckpointV3AuthorityGradient>
decode_rows(runtime::checkpoint_v2::Decoder &decoder) {
  const std::uint64_t count = decoder.u64();
  if (count > kMaximumAuthorityRows) {
    throw runtime::Error("Checkpoint v3 authority row count is invalid");
  }
  std::vector<detail::CheckpointV3AuthorityGradient> result;
  result.reserve(runtime::checkpoint_v2::checked_size(count));
  for (std::uint64_t index = 0U; index < count; ++index) {
    result.push_back({decoder.u64(), decoder.f64()});
  }
  return result;
}

void encode_state_rank(runtime::checkpoint_v2::Encoder &encoder,
                       const detail::CheckpointV3RankAuthority &rank) {
  encoder.i32(rank.rank);
  encoder.i32(rank.owned_box.begin.x);
  encoder.i32(rank.owned_box.begin.y);
  encoder.i32(rank.owned_box.begin.z);
  encoder.i32(rank.owned_box.end.x);
  encoder.i32(rank.owned_box.end.y);
  encoder.i32(rank.owned_box.end.z);
  encoder.string(rank.state_filename);
  encoder.u64(rank.state_logical_bytes);
  encoder.u64(rank.state_actual_bytes);
  encoder.u64(rank.state_crc64);
}

detail::CheckpointV3RankAuthority
decode_state_rank(runtime::checkpoint_v2::Decoder &decoder) {
  detail::CheckpointV3RankAuthority rank;
  rank.rank = decoder.i32();
  rank.owned_box.begin = {decoder.i32(), decoder.i32(), decoder.i32()};
  rank.owned_box.end = {decoder.i32(), decoder.i32(), decoder.i32()};
  rank.state_filename = decoder.string();
  rank.state_logical_bytes = decoder.u64();
  rank.state_actual_bytes = decoder.u64();
  rank.state_crc64 = decoder.u64();
  return rank;
}

std::vector<std::uint8_t> encode_ibm_section(
    const std::vector<detail::CheckpointV3RankAuthority> &ranks) {
  runtime::checkpoint_v2::Encoder encoder;
  encoder.u32(UINT32_C(0x334d4249));
  encoder.u32(static_cast<std::uint32_t>(ranks.size()));
  for (const auto &rank : ranks) {
    encoder.i32(rank.rank);
    for (const auto fingerprint : rank.fingerprints) {
      encoder.u64(fingerprint);
    }
    encoder.boolean(rank.history_available);
    encoder.boolean(rank.committed_available);
    append_rows(encoder, rank.history);
    append_rows(encoder, rank.committed);
  }
  return std::move(encoder).take();
}

void decode_ibm_section(
    const std::vector<std::uint8_t> &bytes,
    std::vector<detail::CheckpointV3RankAuthority> &ranks) {
  runtime::checkpoint_v2::Decoder decoder(bytes);
  if (decoder.u32() != UINT32_C(0x334d4249) ||
      decoder.u32() != static_cast<std::uint32_t>(ranks.size())) {
    throw runtime::Error("Checkpoint v3 IBM section is invalid");
  }
  for (std::size_t index = 0U; index < ranks.size(); ++index) {
    auto &rank = ranks[index];
    if (decoder.i32() != rank.rank) {
      throw runtime::Error("Checkpoint v3 IBM rank identity is invalid");
    }
    for (auto &fingerprint : rank.fingerprints) {
      fingerprint = decoder.u64();
    }
    rank.history_available = decoder.boolean();
    rank.committed_available = decoder.boolean();
    rank.history = decode_rows(decoder);
    rank.committed = decode_rows(decoder);
  }
  decoder.require_eof();
}

std::vector<std::uint8_t> encode_wale_section(
    const detail::CheckpointV3WaleIdentity &identity) {
  runtime::checkpoint_v2::Encoder encoder;
  encoder.u32(UINT32_C(0x334c4157));
  encoder.f64(identity.coefficient);
  encoder.f64(identity.turbulent_prandtl);
  encoder.f64(identity.turbulent_schmidt);
  encoder.u64(identity.numerical_config_crc64);
  encoder.u32(identity.transient_schema_version);
  encoder.u32(static_cast<std::uint32_t>(identity.transient_fields.size()));
  for (const auto field : identity.transient_fields) {
    encoder.u8(static_cast<std::uint8_t>(field));
  }
  return std::move(encoder).take();
}

detail::CheckpointV3WaleIdentity
decode_wale_section(const std::vector<std::uint8_t> &bytes) {
  runtime::checkpoint_v2::Decoder decoder(bytes);
  if (decoder.u32() != UINT32_C(0x334c4157)) {
    throw runtime::Error("Checkpoint v3 WALE section magic is invalid");
  }
  detail::CheckpointV3WaleIdentity result;
  result.coefficient = decoder.f64();
  result.turbulent_prandtl = decoder.f64();
  result.turbulent_schmidt = decoder.f64();
  result.numerical_config_crc64 = decoder.u64();
  result.transient_schema_version = decoder.u32();
  if (decoder.u32() !=
      static_cast<std::uint32_t>(result.transient_fields.size())) {
    throw runtime::Error("Checkpoint v3 WALE field count is invalid");
  }
  for (auto &field : result.transient_fields) {
    field = static_cast<detail::CheckpointV3WaleTransientField>(decoder.u8());
  }
  decoder.require_eof();
  if (!valid_wale_identity(result)) {
    throw runtime::Error("Checkpoint v3 WALE identity is invalid");
  }
  return result;
}

std::vector<std::uint8_t>
encode_ideal_gas_section(const IdealGasClosureState &state) {
  runtime::checkpoint_v2::Encoder encoder;
  encoder.u32(UINT32_C(0x33474949));
  encoder.u8(static_cast<std::uint8_t>(state.mode));
  encoder.f64(state.thermodynamic_pressure_pa);
  encoder.boolean(state.target_mass_kg.has_value());
  encoder.f64(state.target_mass_kg.value_or(0.0));
  encoder.u64(state.revision);
  return std::move(encoder).take();
}

IdealGasClosureState
decode_ideal_gas_section(const std::vector<std::uint8_t> &bytes) {
  runtime::checkpoint_v2::Decoder decoder(bytes);
  if (decoder.u32() != UINT32_C(0x33474949)) {
    throw runtime::Error("Checkpoint v3 ideal-gas section magic is invalid");
  }
  IdealGasClosureState result;
  result.mode = static_cast<IdealGasPressureMode>(decoder.u8());
  result.thermodynamic_pressure_pa = decoder.f64();
  const bool has_target = decoder.boolean();
  const double target = decoder.f64();
  if (has_target) {
    result.target_mass_kg = target;
  } else if (bits(target) != bits(0.0)) {
    throw runtime::Error(
        "Checkpoint v3 absent ideal-gas target is noncanonical");
  }
  result.revision = decoder.u64();
  decoder.require_eof();
  if (!valid_ideal_gas_state(result)) {
    throw runtime::Error("Checkpoint v3 ideal-gas state is invalid");
  }
  return result;
}

std::vector<std::uint8_t>
encode_rank_authority(const detail::CheckpointV3RankAuthority &rank) {
  if (!valid_rank_authority(rank)) {
    throw runtime::Error("Checkpoint v3 rank authority is invalid");
  }
  runtime::checkpoint_v2::Encoder encoder;
  encoder.u32(UINT32_C(0x33415256));
  encoder.i32(rank.rank);
  encoder.i32(rank.owned_box.begin.x);
  encoder.i32(rank.owned_box.begin.y);
  encoder.i32(rank.owned_box.begin.z);
  encoder.i32(rank.owned_box.end.x);
  encoder.i32(rank.owned_box.end.y);
  encoder.i32(rank.owned_box.end.z);
  for (const auto fingerprint : rank.fingerprints) {
    encoder.u64(fingerprint);
  }
  encoder.string(rank.state_filename);
  encoder.u64(rank.state_logical_bytes);
  encoder.u64(rank.state_actual_bytes);
  encoder.u64(rank.state_crc64);
  encoder.boolean(rank.history_available);
  encoder.boolean(rank.committed_available);
  append_rows(encoder, rank.history);
  append_rows(encoder, rank.committed);
  return std::move(encoder).take();
}

detail::CheckpointV3RankAuthority
decode_rank_authority(const std::vector<std::uint8_t> &bytes) {
  runtime::checkpoint_v2::Decoder decoder(bytes);
  if (decoder.u32() != UINT32_C(0x33415256)) {
    throw runtime::Error("Checkpoint v3 rank authority magic is invalid");
  }
  detail::CheckpointV3RankAuthority rank;
  rank.rank = decoder.i32();
  rank.owned_box.begin = {decoder.i32(), decoder.i32(), decoder.i32()};
  rank.owned_box.end = {decoder.i32(), decoder.i32(), decoder.i32()};
  for (auto &fingerprint : rank.fingerprints) {
    fingerprint = decoder.u64();
  }
  rank.state_filename = decoder.string();
  rank.state_logical_bytes = decoder.u64();
  rank.state_actual_bytes = decoder.u64();
  rank.state_crc64 = decoder.u64();
  rank.history_available = decoder.boolean();
  rank.committed_available = decoder.boolean();
  rank.history = decode_rows(decoder);
  rank.committed = decode_rows(decoder);
  decoder.require_eof();
  if (!valid_rank_authority(rank)) {
    throw runtime::Error("Checkpoint v3 rank authority is invalid");
  }
  return rank;
}

std::vector<std::vector<std::uint8_t>> allgather_variable_bytes(
    const runtime::MpiContext &mpi, const std::vector<std::uint8_t> &local,
    std::uint64_t &collectives) {
  auto size_status = runtime::checkpoint_v2::converge_phase(
      mpi, local.size() <= static_cast<std::size_t>(INT_MAX), collectives,
      "MPI_Allreduce(Checkpoint v3 authority sizes)");
  if (!size_status.ok) {
    throw runtime::checkpoint_v2::CollectivePreparationError(
        size_status.failing_rank,
        "Checkpoint v3 authority record exceeds MPI count range");
  }
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts;
  std::vector<int> offsets;
  bool workspace_ok = true;
  try {
    counts.resize(static_cast<std::size_t>(mpi.size()));
    offsets.resize(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    workspace_ok = false;
  }
  const auto workspace_status = runtime::checkpoint_v2::converge_phase(
      mpi, workspace_ok, collectives,
      "MPI_Allreduce(Checkpoint v3 authority workspace)");
  if (!workspace_status.ok) {
    throw runtime::checkpoint_v2::CollectivePreparationError(
        workspace_status.failing_rank,
        "Checkpoint v3 authority workspace is unavailable");
  }
  if (MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    mpi.comm()) != MPI_SUCCESS) {
    throw runtime::Error("MPI_Allgather(Checkpoint v3 authority counts) failed");
  }
  ++collectives;
  std::size_t total{};
  bool layout_ok = true;
  for (std::size_t index = 0U; index < counts.size(); ++index) {
    if (counts[index] < 0 || total > static_cast<std::size_t>(INT_MAX) ||
        static_cast<std::size_t>(counts[index]) >
            static_cast<std::size_t>(INT_MAX) - total) {
      layout_ok = false;
      break;
    }
    offsets[index] = static_cast<int>(total);
    total += static_cast<std::size_t>(counts[index]);
  }
  std::vector<std::uint8_t> gathered;
  bool allocation_ok = layout_ok;
  if (allocation_ok) {
    try {
      gathered.resize(total);
    } catch (...) {
      allocation_ok = false;
    }
  }
  const auto allocation_status = runtime::checkpoint_v2::converge_phase(
      mpi, allocation_ok, collectives,
      "MPI_Allreduce(Checkpoint v3 authority buffer)");
  if (!allocation_status.ok) {
    throw runtime::checkpoint_v2::CollectivePreparationError(
        allocation_status.failing_rank,
        "Checkpoint v3 authority gather buffer is unavailable");
  }
  if (MPI_Allgatherv(local.data(), local_count, MPI_BYTE, gathered.data(),
                     counts.data(), offsets.data(), MPI_BYTE, mpi.comm()) !=
      MPI_SUCCESS) {
    throw runtime::Error("MPI_Allgatherv(Checkpoint v3 authority) failed");
  }
  ++collectives;
  std::vector<std::vector<std::uint8_t>> result;
  bool result_ok = true;
  try {
    result.reserve(counts.size());
    for (std::size_t index = 0U; index < counts.size(); ++index) {
      const auto begin = gathered.begin() + offsets[index];
      result.emplace_back(begin, begin + counts[index]);
    }
  } catch (...) {
    result_ok = false;
  }
  const auto result_status = runtime::checkpoint_v2::converge_phase(
      mpi, result_ok, collectives,
      "MPI_Allreduce(Checkpoint v3 authority records)");
  if (!result_status.ok) {
    throw runtime::checkpoint_v2::CollectivePreparationError(
        result_status.failing_rank,
        "Checkpoint v3 authority records are unavailable");
  }
  return result;
}

bool same_rows(
    const std::vector<detail::CheckpointV3AuthorityGradient> &left,
    const std::vector<detail::CheckpointV3AuthorityGradient> &right) noexcept {
  if (left.size() != right.size()) {
    return false;
  }
  for (std::size_t index = 0U; index < left.size(); ++index) {
    if (left[index].link != right[index].link ||
        bits(left[index].value) != bits(right[index].value)) {
      return false;
    }
  }
  return true;
}

} // namespace

struct detail::CheckpointV3Access final {
  static CheckpointV3Report make(
      CheckpointV3Operation operation, CheckpointV3Disposition disposition,
      CheckpointV3FailureReason reason, CheckpointV3Phase phase, int rank,
      std::uint64_t step, double time_s, std::uint64_t manifest_crc,
      CheckpointV3CheckStatus crc, CheckpointV3CheckStatus fingerprint,
      CheckpointV3CheckStatus partition,
      CheckpointV3CheckStatus rollback) noexcept {
    CheckpointV3Report result;
    result.operation_ = operation;
    result.disposition_ = disposition;
    result.reason_ = reason;
    result.phase_ = phase;
    result.lowest_failing_rank_ = rank;
    result.step_ = step;
    result.time_s_ = time_s;
    result.manifest_crc64_ = manifest_crc;
    result.crc_status_ = crc;
    result.fingerprint_status_ = fingerprint;
    result.partition_status_ = partition;
    result.rollback_status_ = rollback;
    return result;
  }
  static CheckpointV3ReadResult make(CheckpointV3Report report,
                                     CheckpointV3ControlState control,
                                     bool restored) noexcept {
    CheckpointV3ReadResult result;
    result.report_ = std::move(report);
    result.control_ = control;
    result.restored_ = restored;
    return result;
  }
  static CheckpointV3Report
  with_presence(CheckpointV3Report report,
                CheckpointV3Presence presence) noexcept {
    report.presence_ = presence;
    return report;
  }
  static CheckpointV3ReadResult
  with_presence(CheckpointV3ReadResult result,
                CheckpointV3Presence presence) noexcept {
    result.report_.presence_ = presence;
    return result;
  }
};

namespace {

CheckpointV3Report failed_report(CheckpointV3Operation operation,
                                 CheckpointV3FailureReason reason,
                                 CheckpointV3Phase phase, int rank,
                                 const FlowState &state,
                                 CheckpointV3CheckStatus crc =
                                     CheckpointV3CheckStatus::not_checked,
                                 CheckpointV3CheckStatus fingerprint =
                                     CheckpointV3CheckStatus::not_checked,
                                 CheckpointV3CheckStatus partition =
                                     CheckpointV3CheckStatus::not_checked,
                                 CheckpointV3CheckStatus rollback =
                                     CheckpointV3CheckStatus::not_checked) {
  AcceptedStepMetadata metadata{};
  if (detail::FlowStateCheckpointAccess::live(state)) {
    metadata = state.metadata();
  }
  return detail::CheckpointV3Access::make(
      operation, CheckpointV3Disposition::failed, reason, phase, rank,
      metadata.step, metadata.time_s, 0U, crc, fingerprint, partition,
      rollback);
}

CheckpointV3ReadResult failed_read_result(
    CheckpointV3FailureReason reason, CheckpointV3Phase phase, int rank,
    const FlowState &state,
    CheckpointV3CheckStatus crc = CheckpointV3CheckStatus::not_checked,
    CheckpointV3CheckStatus fingerprint =
        CheckpointV3CheckStatus::not_checked,
    CheckpointV3CheckStatus partition =
        CheckpointV3CheckStatus::not_checked) {
  return detail::CheckpointV3Access::make(
      failed_report(CheckpointV3Operation::read, reason, phase, rank, state,
                    crc, fingerprint, partition,
                    CheckpointV3CheckStatus::passed),
      {}, false);
}

std::uint64_t checkpoint_numerical_config_crc64(
    const config::ImmersedFlowCaseConfig &config) {
  auto numerical_config = config;
  numerical_config.common_flow.case_name = "checkpoint-v3-numerical-identity";
  numerical_config.common_flow.resources.expected_ranks.reset();
  numerical_config.common_flow.resources.process_grid.reset();
  numerical_config.common_flow.time.steps = 1;
  numerical_config.common_flow.restart.read = false;
  numerical_config.common_flow.restart.read_directory.reset();
  numerical_config.common_flow.restart.write_directory = "checkpoint";
  numerical_config.common_flow.restart.write_interval = 1;
  numerical_config.common_flow.diagnostics.directory = "diagnostics";
  numerical_config.common_flow.diagnostics.write_interval = 1;
  numerical_config.common_flow.diagnostics.write_mesh = false;
  numerical_config.common_flow.performance.enabled = false;
  numerical_config.common_flow.performance.directory = "performance";
  numerical_config.common_flow.performance.warmup_steps = 1;
  numerical_config.common_flow.performance.measured_steps = 1;
  numerical_config.common_flow.performance.repetitions = 1;
  if (numerical_config.immersed_boundary.geometry) {
    numerical_config.immersed_boundary.geometry->file = "surface.stl";
  }
  const auto canonical = config::to_resolved_json_v3(
      config::ResolvedCaseV3(std::move(numerical_config)));
  return runtime::checkpoint_v2::crc64_ecma(canonical.data(),
                                             canonical.size());
}

std::array<std::uint64_t, 9> checkpoint_identity(
    const config::ImmersedFlowCaseConfig &config,
    const immersed::ImmersedSurface &surface,
    const immersed::SurfaceQuery &query, const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost,
    const immersed::WallQuadraturePlan &wall,
    const immersed::LocalFlowPatternTransform &transform) {
  return {surface.fingerprint(),
          query.fingerprint(),
          domain.classification_fingerprint(),
          domain.surface_coverage_fingerprint(),
          domain.active_cells().fingerprint(),
          domain.active_boundaries().fingerprint(),
          ghost.fingerprint(),
          wall.fingerprint(),
          transform.algorithm_fingerprint() ^
              checkpoint_numerical_config_crc64(config)};
}

std::string rank_state_filename(int rank) {
  std::array<char, 32> text{};
  const int count =
      std::snprintf(text.data(), text.size(), "rank-%06d.v3.bin", rank);
  if (count <= 0 || static_cast<std::size_t>(count) >= text.size()) {
    throw runtime::Error("Checkpoint v3 rank filename is invalid");
  }
  return std::string(text.data(), static_cast<std::size_t>(count));
}

std::vector<std::string> checkpoint_inventory(int ranks) {
  std::vector<std::string> result{"COMPLETED", "manifest.v3.bin"};
  for (int rank = 0; rank < ranks; ++rank) {
    result.push_back(rank_state_filename(rank));
  }
  return result;
}

bool exact_checkpoint_inventory(const std::filesystem::path &directory,
                                int ranks) {
  return runtime::checkpoint_v2::exact_directory_inventory(
      directory, checkpoint_inventory(ranks));
}

bool finite_values(const std::vector<double> &values) noexcept {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

bool valid_checkpoint_layer(const FlowLayerValues &layer,
                            const mesh::MeshTopology &topology,
                            const immersed::ImmersedDomain *domain,
                            config::DensityModel density_model) noexcept {
  const bool constant = density_model == config::DensityModel::constant;
  const bool material = density_model == config::DensityModel::material;
  const bool ideal_gas = density_model == config::DensityModel::ideal_gas;
  if (layer.density.size() != topology.owned_cell_count() ||
      layer.velocity.size() != topology.owned_cell_count() * 3U ||
      layer.mechanical_pressure.size() != topology.owned_cell_count() ||
      layer.face_velocity.size() != topology.local_face_count() * 3U ||
      layer.face_mass_flux.size() != topology.local_face_count() ||
      (!constant && !material && !ideal_gas) ||
      (constant && !layer.transported_cell_fields.empty()) ||
      (!constant && layer.transported_cell_fields.empty()) ||
      !finite_values(layer.density) || !finite_values(layer.velocity) ||
      !finite_values(layer.mechanical_pressure) ||
      !finite_values(layer.face_velocity) ||
      !finite_values(layer.face_mass_flux)) {
    return false;
  }
  if (!std::all_of(layer.transported_cell_fields.begin(),
                   layer.transported_cell_fields.end(),
                   [&](const std::vector<double> &values) {
                     return values.size() == topology.owned_cell_count() &&
                            finite_values(values);
                   })) {
    return false;
  }
  const auto positive_zero = [](double value) noexcept {
    return bits(value) == bits(0.0);
  };
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const bool fluid = domain == nullptr ||
                       domain->region(cell) == immersed::CellRegion::fluid;
    if (fluid) {
      if (!(layer.density[cell] > 0.0) ||
          (ideal_gas && !(layer.transported_cell_fields.front()[cell] > 0.0))) {
        return false;
      }
      continue;
    }
    if (!positive_zero(layer.density[cell]) ||
        !positive_zero(layer.mechanical_pressure[cell])) {
      return false;
    }
    for (std::size_t component = 0U; component < 3U; ++component) {
      if (!positive_zero(layer.velocity[cell * 3U + component])) {
        return false;
      }
    }
    for (const auto &transported : layer.transported_cell_fields) {
      if (!positive_zero(transported[cell])) {
        return false;
      }
    }
  }
  if (domain == nullptr) {
    return true;
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    const bool owner_active =
        domain->region(topology.owner(face)) == immersed::CellRegion::fluid;
    const bool neighbour_active =
        neighbour.has_value() &&
        domain->region(*neighbour) == immersed::CellRegion::fluid;
    if (owner_active && (!neighbour.has_value() || neighbour_active)) {
      continue;
    }
    if (!positive_zero(layer.face_mass_flux[face])) {
      return false;
    }
    for (std::size_t component = 0U; component < 3U; ++component) {
      if (!positive_zero(layer.face_velocity[face * 3U + component])) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

std::vector<std::uint8_t> detail::encode_checkpoint_v3_manifest(
    const CheckpointV3Manifest &value) {
  require_valid(value);
  runtime::checkpoint_v2::Encoder encoder;
  encoder.raw(kManifestMagic.data(), kManifestMagic.size());
  encoder.u32(kManifestSchema);
  encoder.u8(static_cast<std::uint8_t>(value.presence));
  encoder.i32(value.rank_count);
  encoder.i32(value.process_grid.x);
  encoder.i32(value.process_grid.y);
  encoder.i32(value.process_grid.z);
  encoder.u64(value.payload_report_fingerprint);
  encoder.u64(value.payload_manifest_crc64);
  encoder.u64(value.metadata.step);
  encoder.f64(value.metadata.time_s);
  encoder.f64(value.metadata.dt_s);
  encoder.f64(value.metadata.previous_dt_s);
  encoder.u8(static_cast<std::uint8_t>(value.metadata.order));
  encoder.f64(value.control.proposed_next_dt_s);
  encoder.u32(value.control.last_retry_count);
  encoder.u32(static_cast<std::uint32_t>(value.ranks.size()));
  if (value.presence == CheckpointV3Presence::constant_static_ibm) {
    for (const auto &rank : value.ranks) {
      encoder.i32(rank.rank);
      encoder.i32(rank.owned_box.begin.x);
      encoder.i32(rank.owned_box.begin.y);
      encoder.i32(rank.owned_box.begin.z);
      encoder.i32(rank.owned_box.end.x);
      encoder.i32(rank.owned_box.end.y);
      encoder.i32(rank.owned_box.end.z);
      for (const auto fingerprint : rank.fingerprints) {
        encoder.u64(fingerprint);
      }
      encoder.string(rank.state_filename);
      encoder.u64(rank.state_logical_bytes);
      encoder.u64(rank.state_actual_bytes);
      encoder.u64(rank.state_crc64);
      encoder.boolean(rank.history_available);
      encoder.boolean(rank.committed_available);
      append_rows(encoder, rank.history);
      append_rows(encoder, rank.committed);
    }
  } else {
    for (const auto &rank : value.ranks) {
      encode_state_rank(encoder, rank);
    }
    encoder.u32(value.ibm_section_count);
    encoder.u64(value.ibm_section_bytes);
    if (value.ibm_section_count == 1U) {
      const auto section = encode_ibm_section(value.ranks);
      encoder.raw(section.data(), section.size());
    }
    encoder.u32(value.wale_section_count);
    encoder.u64(value.wale_section_bytes);
    if (value.wale_section_count == 1U) {
      const auto section = encode_wale_section(*value.wale);
      encoder.raw(section.data(), section.size());
    }
    if (profile_has_material_fields(value.presence)) {
      encoder.u32(value.ideal_gas_section_count);
      encoder.u64(value.ideal_gas_section_bytes);
      if (value.ideal_gas_section_count == 1U) {
        const auto section = encode_ideal_gas_section(*value.ideal_gas);
        encoder.raw(section.data(), section.size());
      }
    }
  }
  return std::move(encoder).take();
}

detail::CheckpointV3Manifest detail::decode_checkpoint_v3_manifest(
    const std::vector<std::uint8_t> &bytes) {
  runtime::checkpoint_v2::Decoder decoder(bytes);
  if (decoder.raw(kManifestMagic.size()) !=
      std::vector<std::uint8_t>(kManifestMagic.begin(), kManifestMagic.end())) {
    throw runtime::Error("Checkpoint v3 manifest magic is invalid");
  }
  if (decoder.u32() != kManifestSchema) {
    throw runtime::Error("Checkpoint v3 manifest schema is invalid");
  }
  CheckpointV3Manifest result;
  result.presence = static_cast<CheckpointV3Presence>(decoder.u8());
  if (!known_presence(result.presence)) {
    throw runtime::Error("Checkpoint v3 manifest presence is invalid");
  }
  result.rank_count = decoder.i32();
  result.process_grid = {decoder.i32(), decoder.i32(), decoder.i32()};
  result.payload_report_fingerprint = decoder.u64();
  result.payload_manifest_crc64 = decoder.u64();
  result.metadata.step = decoder.u64();
  result.metadata.time_s = decoder.f64();
  result.metadata.dt_s = decoder.f64();
  result.metadata.previous_dt_s = decoder.f64();
  result.metadata.order = static_cast<MomentumTimeOrder>(decoder.u8());
  result.control.proposed_next_dt_s = decoder.f64();
  result.control.last_retry_count = decoder.u32();
  const std::uint32_t authority_ranks = decoder.u32();
  if (authority_ranks > static_cast<std::uint32_t>(
                            std::numeric_limits<std::int32_t>::max())) {
    throw runtime::Error("Checkpoint v3 authority rank count is invalid");
  }
  result.ranks.reserve(authority_ranks);
  if (result.presence == CheckpointV3Presence::constant_static_ibm) {
    for (std::uint32_t index = 0U; index < authority_ranks; ++index) {
      CheckpointV3RankAuthority rank;
      rank.rank = decoder.i32();
      rank.owned_box.begin = {decoder.i32(), decoder.i32(), decoder.i32()};
      rank.owned_box.end = {decoder.i32(), decoder.i32(), decoder.i32()};
      for (auto &fingerprint : rank.fingerprints) {
        fingerprint = decoder.u64();
      }
      rank.state_filename = decoder.string();
      rank.state_logical_bytes = decoder.u64();
      rank.state_actual_bytes = decoder.u64();
      rank.state_crc64 = decoder.u64();
      rank.history_available = decoder.boolean();
      rank.committed_available = decoder.boolean();
      rank.history = decode_rows(decoder);
      rank.committed = decode_rows(decoder);
      result.ranks.push_back(std::move(rank));
    }
  } else {
    for (std::uint32_t index = 0U; index < authority_ranks; ++index) {
      result.ranks.push_back(decode_state_rank(decoder));
    }
    result.ibm_section_count = decoder.u32();
    result.ibm_section_bytes = decoder.u64();
    if (result.ibm_section_count > 1U ||
        (result.ibm_section_count == 0U &&
         result.ibm_section_bytes != 0U)) {
      throw runtime::Error("Checkpoint v3 IBM section presence is invalid");
    }
    if (result.ibm_section_count == 1U) {
      auto section = decoder.raw(
          runtime::checkpoint_v2::checked_size(result.ibm_section_bytes));
      decode_ibm_section(section, result.ranks);
    }
    result.wale_section_count = decoder.u32();
    result.wale_section_bytes = decoder.u64();
    if (result.wale_section_count > 1U ||
        (result.wale_section_count == 0U &&
         result.wale_section_bytes != 0U) ||
        (result.wale_section_count == 1U &&
         result.wale_section_bytes != kWaleSectionBytes)) {
      throw runtime::Error("Checkpoint v3 WALE section presence is invalid");
    }
    if (result.wale_section_count == 1U) {
      result.wale = decode_wale_section(decoder.raw(
          runtime::checkpoint_v2::checked_size(result.wale_section_bytes)));
    }
    if (profile_has_material_fields(result.presence)) {
      result.ideal_gas_section_count = decoder.u32();
      result.ideal_gas_section_bytes = decoder.u64();
      if (result.ideal_gas_section_count > 1U ||
          (result.ideal_gas_section_count == 0U &&
           result.ideal_gas_section_bytes != 0U) ||
          (result.ideal_gas_section_count == 1U &&
           result.ideal_gas_section_bytes != kIdealGasSectionBytes)) {
        throw runtime::Error(
            "Checkpoint v3 ideal-gas section presence is invalid");
      }
      if (result.ideal_gas_section_count == 1U) {
        result.ideal_gas = decode_ideal_gas_section(decoder.raw(
            runtime::checkpoint_v2::checked_size(
                result.ideal_gas_section_bytes)));
      }
    }
  }
  decoder.require_eof();
  require_valid(result);
  return result;
}

bool detail::checkpoint_v3_manifest_equal(
    const CheckpointV3Manifest &left,
    const CheckpointV3Manifest &right) noexcept {
  const auto same_wale = [](const auto &a, const auto &b) {
    if (a.has_value() != b.has_value()) {
      return false;
    }
    return !a.has_value() ||
           (bits(a->coefficient) == bits(b->coefficient) &&
            bits(a->turbulent_prandtl) == bits(b->turbulent_prandtl) &&
            bits(a->turbulent_schmidt) == bits(b->turbulent_schmidt) &&
            a->numerical_config_crc64 == b->numerical_config_crc64 &&
            a->transient_schema_version == b->transient_schema_version &&
            a->transient_fields == b->transient_fields);
  };
  const auto same_ideal_gas = [](const auto &a, const auto &b) {
    if (a.has_value() != b.has_value()) {
      return false;
    }
    if (!a.has_value()) {
      return true;
    }
    if (a->mode != b->mode ||
        bits(a->thermodynamic_pressure_pa) !=
            bits(b->thermodynamic_pressure_pa) ||
        a->target_mass_kg.has_value() != b->target_mass_kg.has_value() ||
        a->revision != b->revision) {
      return false;
    }
    return !a->target_mass_kg.has_value() ||
           bits(*a->target_mass_kg) == bits(*b->target_mass_kg);
  };
  return left.presence == right.presence &&
         left.rank_count == right.rank_count &&
         same(left.process_grid, right.process_grid) &&
         left.payload_report_fingerprint == right.payload_report_fingerprint &&
         left.payload_manifest_crc64 == right.payload_manifest_crc64 &&
         left.metadata.step == right.metadata.step &&
         bits(left.metadata.time_s) == bits(right.metadata.time_s) &&
         bits(left.metadata.dt_s) == bits(right.metadata.dt_s) &&
         bits(left.metadata.previous_dt_s) ==
             bits(right.metadata.previous_dt_s) &&
         left.metadata.order == right.metadata.order &&
         bits(left.control.proposed_next_dt_s) ==
             bits(right.control.proposed_next_dt_s) &&
         left.control.last_retry_count == right.control.last_retry_count &&
         left.ibm_section_count == right.ibm_section_count &&
         left.ibm_section_bytes == right.ibm_section_bytes &&
         left.wale_section_count == right.wale_section_count &&
         left.wale_section_bytes == right.wale_section_bytes &&
         same_wale(left.wale, right.wale) &&
         left.ideal_gas_section_count == right.ideal_gas_section_count &&
         left.ideal_gas_section_bytes == right.ideal_gas_section_bytes &&
         same_ideal_gas(left.ideal_gas, right.ideal_gas) &&
         left.ranks.size() == right.ranks.size() &&
         std::equal(left.ranks.begin(), left.ranks.end(), right.ranks.begin(),
                    [](const CheckpointV3RankAuthority &a,
                       const CheckpointV3RankAuthority &b) {
                      return a.rank == b.rank &&
                             same(a.owned_box, b.owned_box) &&
                             a.fingerprints == b.fingerprints &&
                             a.state_filename == b.state_filename &&
                             a.state_logical_bytes == b.state_logical_bytes &&
                             a.state_actual_bytes == b.state_actual_bytes &&
                             a.state_crc64 == b.state_crc64 &&
                             a.history_available == b.history_available &&
                             a.committed_available == b.committed_available &&
                             same_rows(a.history, b.history) &&
                             same_rows(a.committed, b.committed);
                    });
}

std::uint64_t detail::checkpoint_v3_manifest_crc64(
    const std::vector<std::uint8_t> &bytes) noexcept {
  return runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size());
}

std::size_t detail::checkpoint_v3_manifest_presence_offset() noexcept {
  return kManifestMagic.size() + sizeof(std::uint32_t);
}

namespace {

bool valid_write_modules(const CheckpointV3WriteModules &modules) noexcept {
  const bool has_ibm = modules.surface != nullptr && modules.query != nullptr &&
                       modules.domain != nullptr &&
                       modules.ghost_plan != nullptr &&
                       modules.wall_plan != nullptr &&
                       modules.transform != nullptr;
  const bool has_no_ibm = modules.surface == nullptr &&
                          modules.query == nullptr &&
                          modules.domain == nullptr &&
                          modules.ghost_plan == nullptr &&
                          modules.wall_plan == nullptr &&
                          modules.transform == nullptr;
  if (!known_presence(modules.presence)) {
    return false;
  }
  const bool expected_ibm = profile_has_ibm(modules.presence);
  return (expected_ibm ? has_ibm : has_no_ibm) &&
         (modules.flow != nullptr) == expected_ibm &&
         (modules.wale != nullptr) == profile_has_wale(modules.presence) &&
         (modules.ideal_gas != nullptr) ==
             profile_has_ideal_gas(modules.presence);
}

bool valid_read_modules(const CheckpointV3ReadModules &modules) noexcept {
  const bool has_ibm = modules.surface != nullptr && modules.query != nullptr &&
                       modules.domain != nullptr &&
                       modules.ghost_plan != nullptr &&
                       modules.wall_plan != nullptr &&
                       modules.transform != nullptr;
  const bool has_no_ibm = modules.surface == nullptr &&
                          modules.query == nullptr &&
                          modules.domain == nullptr &&
                          modules.ghost_plan == nullptr &&
                          modules.wall_plan == nullptr &&
                          modules.transform == nullptr;
  if (!known_presence(modules.presence)) {
    return false;
  }
  const bool expected_ibm = profile_has_ibm(modules.presence);
  return (expected_ibm ? has_ibm : has_no_ibm) &&
         (modules.flow != nullptr) == expected_ibm &&
         (modules.wale != nullptr) == profile_has_wale(modules.presence) &&
         (modules.ideal_gas != nullptr) ==
             profile_has_ideal_gas(modules.presence);
}

detail::CheckpointV3WaleIdentity
wale_identity(const les::WaleModel &wale,
              const config::ImmersedFlowCaseConfig &config) {
  const auto control = wale.control();
  return {control.coefficient,
          control.turbulent_prandtl,
          control.turbulent_schmidt,
          checkpoint_numerical_config_crc64(config),
          kWaleTransientSchemaVersion,
          {detail::CheckpointV3WaleTransientField::nu_t_m2_per_s,
           detail::CheckpointV3WaleTransientField::mu_sgs_pa_s,
           detail::CheckpointV3WaleTransientField::mu_eff_pa_s}};
}

bool same_wale_identity(const detail::CheckpointV3WaleIdentity &left,
                        const detail::CheckpointV3WaleIdentity &right) noexcept {
  return bits(left.coefficient) == bits(right.coefficient) &&
         bits(left.turbulent_prandtl) == bits(right.turbulent_prandtl) &&
         bits(left.turbulent_schmidt) == bits(right.turbulent_schmidt) &&
         left.numerical_config_crc64 == right.numerical_config_crc64 &&
         left.transient_schema_version == right.transient_schema_version &&
         left.transient_fields == right.transient_fields;
}

config::DensityModel
profile_density_model(CheckpointV3Presence presence) noexcept {
  switch (presence) {
  case CheckpointV3Presence::constant_static_ibm:
  case CheckpointV3Presence::constant_body_fitted_wale:
  case CheckpointV3Presence::constant_static_ibm_wale:
    return config::DensityModel::constant;
  case CheckpointV3Presence::material_static_ibm:
  case CheckpointV3Presence::material_body_fitted_wale:
  case CheckpointV3Presence::material_static_ibm_wale:
    return config::DensityModel::material;
  case CheckpointV3Presence::ideal_gas_static_ibm:
  case CheckpointV3Presence::ideal_gas_body_fitted_wale:
  case CheckpointV3Presence::ideal_gas_static_ibm_wale:
    return config::DensityModel::ideal_gas;
  }
  return config::DensityModel::constant;
}

bool supported_profile(
    const runtime::MpiContext &mpi,
    const config::ImmersedFlowCaseConfig &config,
    CheckpointV3Presence presence, const les::WaleModel *wale) noexcept {
  const bool ibm = profile_has_ibm(presence);
  const bool has_wale = profile_has_wale(presence);
  if (mpi.size() <= 0 || config.schema_version != 3 ||
      config.common_flow.schema_version != 2 ||
      config.common_flow.simulation_type !=
          config::SimulationType::variable_density_flow ||
      config.common_flow.density_model != profile_density_model(presence) ||
      config.common_flow.time.mode != config::TimeMode::fixed ||
      (profile_density_model(presence) == config::DensityModel::constant &&
       !config.common_flow.scalars.empty()) ||
      (has_wale && (config.les.model != config::LesModel::wale ||
                    !config.les.wale.has_value() || wale == nullptr)) ||
      (!has_wale && (config.les.model != config::LesModel::none ||
                     config.les.wale.has_value() || wale != nullptr)) ||
      (ibm &&
       (config.immersed_boundary.model !=
            config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell ||
        !config.immersed_boundary.geometry.has_value() ||
        !config.immersed_boundary.wall.has_value())) ||
      (!ibm &&
       (config.immersed_boundary.model != config::ImmersedBoundaryModel::none ||
        config.immersed_boundary.geometry.has_value() ||
        config.immersed_boundary.wall.has_value()))) {
    return false;
  }
  if (!has_wale) {
    return true;
  }
  const auto control = wale->control();
  return bits(control.coefficient) == bits(config.les.wale->coefficient) &&
         bits(control.turbulent_prandtl) ==
             bits(config.les.wale->turbulent_prandtl) &&
         bits(control.turbulent_schmidt) ==
             bits(config.les.wale->turbulent_schmidt);
}

bool supported_write_modules(
    const runtime::MpiContext &mpi,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3WriteModules &modules,
    const CheckpointV3ControlState &control) noexcept {
  if (!valid_write_modules(modules)) {
    return false;
  }
  return supported_profile(mpi, config, modules.presence, modules.wale) &&
         std::isfinite(control.proposed_next_dt_s) &&
         control.proposed_next_dt_s == config.common_flow.time.initial_dt_s &&
         control.last_retry_count <= 8U;
}

bool supported_read_modules(
    const runtime::MpiContext &mpi,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3ReadModules &modules) noexcept {
  if (!valid_read_modules(modules)) {
    return false;
  }
  return supported_profile(mpi, config, modules.presence, modules.wale);
}

CheckpointV3Report with_presence(CheckpointV3Report report,
                                 CheckpointV3Presence presence) noexcept {
  return detail::CheckpointV3Access::with_presence(std::move(report), presence);
}

CheckpointV3ReadResult with_presence(CheckpointV3ReadResult result,
                                     CheckpointV3Presence presence) noexcept {
  return detail::CheckpointV3Access::with_presence(std::move(result), presence);
}

CheckpointV3Report write_checkpoint_v3_modules(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &boundaries,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3WriteModules &modules, const FlowState &state,
    CheckpointV3ControlState control,
    const std::filesystem::path &directory) {
  std::uint64_t collectives{};
  const auto support_status = runtime::checkpoint_v2::converge_phase(
      mpi, supported_write_modules(mpi, config, modules, control), collectives,
      "MPI_Allreduce(Checkpoint v3 write presence)");
  if (!support_status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::presence,
                         CheckpointV3Phase::preflight,
                         support_status.failing_rank, state);
  }
  const auto state_status = runtime::checkpoint_v2::converge_phase(
      mpi, detail::FlowStateCheckpointAccess::live(state) &&
               !detail::FlowStateCheckpointAccess::attempt_active(state),
      collectives, "MPI_Allreduce(Checkpoint v3 write state readiness)");
  if (!state_status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::state,
                         CheckpointV3Phase::preflight,
                         state_status.failing_rank, state);
  }
  std::vector<std::uint8_t> state_bytes;
  std::vector<std::uint8_t> authority_bytes;
  std::vector<std::uint8_t> common_bytes;
  std::optional<IdealGasClosureState> ideal_gas_state;
  detail::CheckpointV3RankAuthority authority;
  bool prepared = true;
  try {
    const auto history = state.snapshot(FlowLayer::history);
    const auto committed = state.snapshot(FlowLayer::committed);
    if (!valid_checkpoint_layer(history, topology, modules.domain,
                                config.common_flow.density_model) ||
        !valid_checkpoint_layer(committed, topology, modules.domain,
                                config.common_flow.density_model)) {
      throw runtime::Error("Checkpoint v3 source state is invalid");
    }
    if (modules.ideal_gas != nullptr) {
      ideal_gas_state =
          detail::IdealGasClosureCheckpointAccess::snapshot(*modules.ideal_gas);
      const bool open = boundaries.open_domain();
      if (!valid_ideal_gas_state(*ideal_gas_state) ||
          ideal_gas_state->mode !=
              (open ? IdealGasPressureMode::open_fixed
                    : IdealGasPressureMode::closed_dynamic) ||
          (open &&
           (!config.common_flow.physics.thermodynamic_pressure_pa.has_value() ||
            bits(ideal_gas_state->thermodynamic_pressure_pa) !=
                bits(*config.common_flow.physics
                          .thermodynamic_pressure_pa)))) {
        throw runtime::Error(
            "Checkpoint v3 ideal-gas source authority is invalid");
      }
    }
    std::uint64_t logical_bytes{};
    state_bytes =
        detail::encode_checkpoint_flow_state_rank_payload(state, logical_bytes);
    record_checkpoint_logical_bytes(logical_bytes);
    const auto state_crc = runtime::checkpoint_v2::crc64_ecma(
        state_bytes.data(), state_bytes.size());
    if (modules.flow != nullptr) {
      authority = detail::ImmersedFlowCheckpointAccess::snapshot(
          *modules.flow, mpi.rank(), decomposition.owned_box());
      authority.fingerprints = checkpoint_identity(
          config, *modules.surface, *modules.query, *modules.domain,
          *modules.ghost_plan, *modules.wall_plan, *modules.transform);
    } else {
      authority.rank = mpi.rank();
      authority.owned_box = decomposition.owned_box();
    }
    authority.state_filename = rank_state_filename(mpi.rank());
    authority.state_logical_bytes = logical_bytes;
    authority.state_actual_bytes = state_bytes.size();
    authority.state_crc64 = state_crc;
    authority_bytes = encode_rank_authority(authority);
    runtime::checkpoint_v2::Encoder common;
    common.string(config::to_resolved_json_v3(config::ResolvedCaseV3(config)));
    common.string(directory.generic_string());
    const auto metadata = state.metadata();
    common.u64(metadata.step);
    common.f64(metadata.time_s);
    common.f64(metadata.dt_s);
    common.f64(metadata.previous_dt_s);
    common.u8(static_cast<std::uint8_t>(metadata.order));
    common.f64(control.proposed_next_dt_s);
    common.u32(control.last_retry_count);
    common.i32(decomposition.process_grid().x);
    common.i32(decomposition.process_grid().y);
    common.i32(decomposition.process_grid().z);
    if (ideal_gas_state.has_value()) {
      common.u8(static_cast<std::uint8_t>(ideal_gas_state->mode));
      common.f64(ideal_gas_state->thermodynamic_pressure_pa);
      common.boolean(ideal_gas_state->target_mass_kg.has_value());
      common.f64(ideal_gas_state->target_mass_kg.value_or(0.0));
      common.u64(ideal_gas_state->revision);
    }
    common_bytes = std::move(common).take();
  } catch (...) {
    prepared = false;
  }
  auto status = runtime::checkpoint_v2::converge_phase(
      mpi, prepared, collectives,
      "MPI_Allreduce(Checkpoint v3 write preparation)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::state,
                         CheckpointV3Phase::rank_payload,
                         status.failing_rank, state);
  }
  try {
    status = runtime::checkpoint_v2::opaque_bytes_agreement(
        mpi, common_bytes, collectives,
        "MPI_Allreduce(Checkpoint v3 common authority)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::state,
                         CheckpointV3Phase::preflight,
                         error.failing_rank(), state);
  }
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::invalid_input,
                         CheckpointV3Phase::preflight,
                         status.failing_rank, state);
  }

  bool directory_ok = true;
  if (mpi.rank() == 0) {
    try {
      runtime::checkpoint_v2::create_directory_exclusive(directory);
    } catch (...) {
      directory_ok = false;
    }
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, directory_ok, collectives,
      "MPI_Allreduce(Checkpoint v3 directory)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::filesystem,
                         CheckpointV3Phase::rank_payload,
                         status.failing_rank, state);
  }

  const auto rank_temp = directory / (authority.state_filename + ".tmp");
  const auto rank_final = directory / authority.state_filename;
  bool rank_written = true;
  try {
    const auto verified = runtime::checkpoint_v2::write_verified_temporary(
        rank_temp, state_bytes);
    rank_written = verified.crc64 == authority.state_crc64 &&
                   verified.actual_size == state_bytes.size();
  } catch (...) {
    rank_written = false;
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, rank_written, collectives,
      "MPI_Allreduce(Checkpoint v3 rank temporary)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::filesystem,
                         CheckpointV3Phase::rank_payload,
                         status.failing_rank, state);
  }
  bool rank_published = true;
  try {
    runtime::checkpoint_v2::publish_no_overwrite(rank_temp, rank_final);
  } catch (...) {
    rank_published = false;
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, rank_published, collectives,
      "MPI_Allreduce(Checkpoint v3 rank publish)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::filesystem,
                         CheckpointV3Phase::rank_payload,
                         status.failing_rank, state);
  }

  std::vector<std::vector<std::uint8_t>> gathered;
  try {
    gathered = allgather_variable_bytes(mpi, authority_bytes, collectives);
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::state,
                         CheckpointV3Phase::manifest,
                         error.failing_rank(), state);
  } catch (...) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::collective_operation,
                         CheckpointV3Phase::manifest, mpi.rank(), state);
  }

  detail::CheckpointV3Manifest manifest;
  bool manifest_prepared = true;
  try {
    manifest.presence = modules.presence;
    manifest.rank_count = mpi.size();
    manifest.process_grid = decomposition.process_grid();
    runtime::checkpoint_v2::Encoder aggregate;
    for (const auto &record : gathered) {
      aggregate.u64(record.size());
      aggregate.raw(record.data(), record.size());
      manifest.ranks.push_back(decode_rank_authority(record));
    }
    const auto aggregate_bytes = std::move(aggregate).take();
    manifest.payload_manifest_crc64 = runtime::checkpoint_v2::crc64_ecma(
        aggregate_bytes.data(), aggregate_bytes.size());
    manifest.payload_report_fingerprint =
        manifest.payload_manifest_crc64 ^ UINT64_C(0x48554e4456335250);
    manifest.metadata = state.metadata();
    manifest.control = control;
    if (profile_has_ibm(modules.presence) &&
        modules.presence != CheckpointV3Presence::constant_static_ibm) {
      const auto bytes = ibm_section_bytes(manifest.ranks);
      if (!bytes.has_value()) {
        throw runtime::Error("Checkpoint v3 IBM section size is invalid");
      }
      manifest.ibm_section_count = 1U;
      manifest.ibm_section_bytes = *bytes;
    }
    if (modules.wale != nullptr) {
      manifest.wale_section_count = 1U;
      manifest.wale_section_bytes = kWaleSectionBytes;
      manifest.wale = wale_identity(*modules.wale, config);
    }
    if (ideal_gas_state.has_value()) {
      manifest.ideal_gas_section_count = 1U;
      manifest.ideal_gas_section_bytes = kIdealGasSectionBytes;
      manifest.ideal_gas = *ideal_gas_state;
    }
    static_cast<void>(detail::encode_checkpoint_v3_manifest(manifest));
  } catch (...) {
    manifest_prepared = false;
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, manifest_prepared, collectives,
      "MPI_Allreduce(Checkpoint v3 manifest preparation)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::state,
                         CheckpointV3Phase::manifest,
                         status.failing_rank, state);
  }

  std::array<std::uint64_t, 2> manifest_authority{};
  bool manifest_written = true;
  if (mpi.rank() == 0) {
    try {
      const auto bytes = detail::encode_checkpoint_v3_manifest(manifest);
      const auto temporary = directory / "manifest.v3.bin.tmp";
      const auto verified = runtime::checkpoint_v2::write_verified_temporary(
          temporary, bytes);
      runtime::checkpoint_v2::publish_no_overwrite(
          temporary, directory / "manifest.v3.bin");
      manifest_authority = {verified.actual_size, verified.crc64};
    } catch (...) {
      manifest_written = false;
    }
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, manifest_written, collectives,
      "MPI_Allreduce(Checkpoint v3 manifest publish)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::filesystem,
                         CheckpointV3Phase::manifest,
                         status.failing_rank, state);
  }
  if (MPI_Bcast(manifest_authority.data(),
                static_cast<int>(manifest_authority.size()), MPI_UINT64_T, 0,
                mpi.comm()) != MPI_SUCCESS) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::collective_operation,
                         CheckpointV3Phase::manifest, mpi.rank(), state);
  }
  ++collectives;

  bool marker_written = true;
  if (mpi.rank() == 0) {
    try {
      const auto bytes = runtime::checkpoint_v2::encode_completed_marker(
          {manifest_authority[0U], manifest_authority[1U],
           manifest.payload_report_fingerprint});
      const auto temporary = directory / "COMPLETED.tmp";
      static_cast<void>(runtime::checkpoint_v2::write_verified_temporary(
          temporary, bytes));
      runtime::checkpoint_v2::publish_no_overwrite(temporary,
                                                    directory / "COMPLETED");
    } catch (...) {
      marker_written = false;
    }
  }
  status = runtime::checkpoint_v2::converge_phase(
      mpi, marker_written, collectives,
      "MPI_Allreduce(Checkpoint v3 completed marker)");
  if (!status.ok) {
    return failed_report(CheckpointV3Operation::write,
                         CheckpointV3FailureReason::filesystem,
                         CheckpointV3Phase::completed_marker,
                         status.failing_rank, state);
  }
  const auto metadata = state.metadata();
  return detail::CheckpointV3Access::make(
      CheckpointV3Operation::write, CheckpointV3Disposition::completed,
      CheckpointV3FailureReason::none, CheckpointV3Phase::completed_marker,
      -1, metadata.step, metadata.time_s, manifest_authority[1U],
      CheckpointV3CheckStatus::passed, CheckpointV3CheckStatus::passed,
      CheckpointV3CheckStatus::passed,
      CheckpointV3CheckStatus::not_checked);
}

CheckpointV3ReadResult read_checkpoint_v3_modules(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3ReadModules &modules, FlowState &state,
    const std::filesystem::path &directory) {
  std::uint64_t collectives{};
  const auto support_status = runtime::checkpoint_v2::converge_phase(
      mpi, supported_read_modules(mpi, config, modules), collectives,
      "MPI_Allreduce(Checkpoint v3 read presence)");
  if (!support_status.ok) {
    return failed_read_result(CheckpointV3FailureReason::presence,
                              CheckpointV3Phase::preflight,
                              support_status.failing_rank, state);
  }

  std::vector<std::uint8_t> path_bytes;
  bool path_prepared = true;
  try {
    const auto normalized = directory.lexically_normal().generic_string();
    if (normalized.empty() || normalized.find('\0') != std::string::npos) {
      throw runtime::Error("Checkpoint v3 read path is invalid");
    }
    runtime::checkpoint_v2::Encoder path_encoder;
    path_encoder.string(normalized);
    path_bytes = std::move(path_encoder).take();
  } catch (...) {
    path_prepared = false;
  }
  const auto path_preparation_status =
      runtime::checkpoint_v2::converge_phase(
          mpi, path_prepared, collectives,
          "MPI_Allreduce(Checkpoint v3 read path preparation)");
  if (!path_preparation_status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::preflight,
                              path_preparation_status.failing_rank, state);
  }
  runtime::checkpoint_v2::CollectiveResult path_status;
  try {
    path_status = runtime::checkpoint_v2::opaque_bytes_agreement(
        mpi, path_bytes, collectives,
        "MPI_Allreduce(Checkpoint v3 read path agreement)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::preflight,
                              error.failing_rank(), state);
  }
  if (!path_status.ok) {
    return failed_read_result(CheckpointV3FailureReason::invalid_input,
                              CheckpointV3Phase::preflight,
                              path_status.failing_rank, state);
  }
  const auto state_status = runtime::checkpoint_v2::converge_phase(
      mpi, detail::FlowStateCheckpointAccess::live(state) &&
               !detail::FlowStateCheckpointAccess::attempt_active(state),
      collectives, "MPI_Allreduce(Checkpoint v3 read state readiness)");
  if (!state_status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::preflight,
                              state_status.failing_rank, state);
  }
  bool closure_ready = true;
  if (modules.ideal_gas != nullptr) {
    try {
      static_cast<void>(
          detail::IdealGasClosureCheckpointAccess::snapshot(*modules.ideal_gas));
    } catch (...) {
      closure_ready = false;
    }
  }
  const auto closure_status = runtime::checkpoint_v2::converge_phase(
      mpi, closure_ready, collectives,
      "MPI_Allreduce(Checkpoint v3 ideal-gas closure readiness)");
  if (!closure_status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::preflight,
                              closure_status.failing_rank, state);
  }

  auto converge = [&](bool local_ok, std::string_view operation) {
    return runtime::checkpoint_v2::converge_phase(mpi, local_ok, collectives,
                                                   operation);
  };
  bool local_ok = true;
  try {
    local_ok = exact_checkpoint_inventory(directory, mpi.size());
  } catch (...) {
    local_ok = false;
  }
  auto status = converge(local_ok,
                         "MPI_Allreduce(Checkpoint v3 read inventory)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::preflight,
                              status.failing_rank, state);
  }

  runtime::checkpoint_v2::CompletedMarker marker;
  std::vector<std::uint8_t> marker_bytes;
  local_ok = true;
  try {
    marker_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "COMPLETED", 40U);
    marker = runtime::checkpoint_v2::decode_completed_marker(marker_bytes);
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok,
                    "MPI_Allreduce(Checkpoint v3 completed marker read)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::completed_marker,
                              status.failing_rank, state);
  }
  try {
    status = runtime::checkpoint_v2::opaque_bytes_agreement(
        mpi, marker_bytes, collectives,
        "MPI_Allreduce(Checkpoint v3 completed marker agreement)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::completed_marker,
                              error.failing_rank(), state);
  }
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::completed_marker,
                              status.failing_rank, state);
  }

  detail::CheckpointV3Manifest manifest;
  std::vector<std::uint8_t> manifest_bytes;
  std::uint64_t manifest_crc{};
  local_ok = true;
  try {
    manifest_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "manifest.v3.bin", marker.manifest_actual_size);
    manifest_crc = runtime::checkpoint_v2::crc64_ecma(
        manifest_bytes.data(), manifest_bytes.size());
    if (manifest_crc != marker.manifest_crc64) {
      throw runtime::Error("Checkpoint v3 manifest CRC is invalid");
    }
    manifest = detail::decode_checkpoint_v3_manifest(manifest_bytes);
    if (marker.common_fingerprint != manifest.payload_report_fingerprint) {
      throw runtime::Error("Checkpoint v3 marker identity is invalid");
    }
    runtime::checkpoint_v2::Encoder aggregate;
    for (const auto &record : manifest.ranks) {
      const auto record_bytes = encode_rank_authority(record);
      aggregate.u64(record_bytes.size());
      aggregate.raw(record_bytes.data(), record_bytes.size());
    }
    const auto aggregate_bytes = std::move(aggregate).take();
    auto aggregate_crc = runtime::checkpoint_v2::crc64_ecma(
        aggregate_bytes.data(), aggregate_bytes.size());
    if (aggregate_crc != manifest.payload_manifest_crc64) {
      throw runtime::Error("Checkpoint v3 rank authority CRC is invalid");
    }
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok, "MPI_Allreduce(Checkpoint v3 manifest read)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::manifest,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::failed);
  }
  try {
    status = runtime::checkpoint_v2::opaque_bytes_agreement(
        mpi, manifest_bytes, collectives,
        "MPI_Allreduce(Checkpoint v3 manifest agreement)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::manifest,
                              error.failing_rank(), state);
  }
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::manifest,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::failed);
  }

  status = converge(manifest.presence == modules.presence,
                    "MPI_Allreduce(Checkpoint v3 profile identity)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::presence,
                              CheckpointV3Phase::manifest,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed);
  }

  const detail::CheckpointV3RankAuthority *rank = nullptr;
  local_ok = true;
  try {
    if (manifest.rank_count != mpi.size() ||
        !same(manifest.process_grid, decomposition.process_grid()) ||
        manifest.ranks.size() != static_cast<std::size_t>(mpi.size())) {
      throw runtime::Error("Checkpoint v3 partition is incompatible");
    }
    rank = &manifest.ranks.at(static_cast<std::size_t>(mpi.rank()));
    if (rank->rank != mpi.rank() ||
        !same(rank->owned_box, decomposition.owned_box()) ||
        rank->state_filename != rank_state_filename(mpi.rank())) {
      throw runtime::Error("Checkpoint v3 rank partition is incompatible");
    }
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok,
                    "MPI_Allreduce(Checkpoint v3 partition identity)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::layout,
                              CheckpointV3Phase::manifest,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::not_checked,
                              CheckpointV3CheckStatus::failed);
  }

  local_ok = true;
  try {
    if (modules.domain != nullptr) {
      const auto expected = checkpoint_identity(
          config, *modules.surface, *modules.query, *modules.domain,
          *modules.ghost_plan, *modules.wall_plan, *modules.transform);
      if (rank->fingerprints != expected) {
        throw runtime::Error("Checkpoint v3 IBM identity is incompatible");
      }
    }
    if (modules.wale != nullptr &&
        (!manifest.wale.has_value() ||
         !same_wale_identity(*manifest.wale,
                             wale_identity(*modules.wale, config)))) {
      throw runtime::Error("Checkpoint v3 WALE identity is incompatible");
    }
    if (manifest.control.proposed_next_dt_s !=
        config.common_flow.time.initial_dt_s) {
      throw runtime::Error("Checkpoint v3 fingerprint is incompatible");
    }
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok,
                    "MPI_Allreduce(Checkpoint v3 geometry identity)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::fingerprint,
                              CheckpointV3Phase::manifest,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::failed,
                              CheckpointV3CheckStatus::passed);
  }

  std::pair<FlowLayerValues, FlowLayerValues> layers;
  local_ok = true;
  try {
    const auto bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / rank->state_filename, rank->state_actual_bytes);
    if (runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size()) !=
        rank->state_crc64) {
      throw runtime::Error("Checkpoint v3 rank CRC is invalid");
    }
    std::uint64_t logical_bytes{};
    layers = detail::decode_checkpoint_flow_state_rank_payload(
        bytes, state, logical_bytes);
    record_checkpoint_logical_bytes(logical_bytes);
    if (logical_bytes != rank->state_logical_bytes ||
        !valid_checkpoint_layer(layers.first, topology, modules.domain,
                                config.common_flow.density_model) ||
        !valid_checkpoint_layer(layers.second, topology, modules.domain,
                                config.common_flow.density_model)) {
      throw runtime::Error("Checkpoint v3 restored state is invalid");
    }
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok, "MPI_Allreduce(Checkpoint v3 rank payload)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::file_integrity,
                              CheckpointV3Phase::rank_payload,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::failed,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed);
  }

  std::optional<FlowState> replacement;
  std::optional<detail::ImmersedFlowCheckpointPreparedRestore>
      authority_replacement;
  std::optional<detail::IdealGasClosureCheckpointPreparedRestore>
      closure_replacement;
  local_ok = true;
  try {
    replacement.emplace(detail::FlowStateCheckpointAccess::prepare_replacement(
        state, layers.first, layers.second, manifest.metadata));
    if (modules.flow != nullptr) {
      authority_replacement.emplace(
          detail::ImmersedFlowCheckpointAccess::prepare_restore(*modules.flow,
                                                                *rank));
    }
    if (modules.ideal_gas != nullptr) {
      if (!manifest.ideal_gas.has_value()) {
        throw runtime::Error(
            "Checkpoint v3 ideal-gas closure authority is absent");
      }
      closure_replacement.emplace(
          detail::IdealGasClosureCheckpointAccess::prepare_restore(
              *modules.ideal_gas, *replacement, *manifest.ideal_gas));
    }
  } catch (...) {
    local_ok = false;
  }
  status = converge(local_ok,
                    "MPI_Allreduce(Checkpoint v3 restore preparation)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::restore_prepare,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed);
  }
  const bool transaction_ready =
      detail::FlowStateCheckpointAccess::read_transaction_ready(state);
  status = converge(transaction_ready,
                    "MPI_Allreduce(Checkpoint v3 restore publish readiness)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::restore_prepare,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed);
  }

  auto report = detail::CheckpointV3Access::make(
      CheckpointV3Operation::read, CheckpointV3Disposition::completed,
      CheckpointV3FailureReason::none, CheckpointV3Phase::restore_publish, -1,
      manifest.metadata.step, manifest.metadata.time_s, manifest_crc,
      CheckpointV3CheckStatus::passed, CheckpointV3CheckStatus::passed,
      CheckpointV3CheckStatus::passed,
      CheckpointV3CheckStatus::not_checked);
  bool entered = true;
  try {
    detail::FlowStateCheckpointAccess::enter_read_transaction(state);
  } catch (...) {
    entered = false;
  }
  status = converge(entered,
                    "MPI_Allreduce(Checkpoint v3 restore transaction entry)");
  if (!status.ok) {
    return failed_read_result(CheckpointV3FailureReason::state,
                              CheckpointV3Phase::restore_prepare,
                              status.failing_rank, state,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed,
                              CheckpointV3CheckStatus::passed);
  }
  detail::FlowStateCheckpointAccess::publish_replacement(
      state, std::move(*replacement));
  if (authority_replacement.has_value()) {
    detail::ImmersedFlowCheckpointAccess::publish_restore(
        *modules.flow, std::move(*authority_replacement), state);
  }
  if (closure_replacement.has_value()) {
    detail::IdealGasClosureCheckpointAccess::publish_restore(
        *modules.ideal_gas, std::move(*closure_replacement));
  }
  return detail::CheckpointV3Access::make(std::move(report), manifest.control,
                                          true);
}

} // namespace

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3WriteModules &modules, const FlowState &state,
    CheckpointV3ControlState control,
    const std::filesystem::path &directory) {
  return with_presence(
      write_checkpoint_v3_modules(mpi, decomposition, topology, geometry,
                                  boundaries, config, modules, state, control,
                                  directory),
      modules.presence);
}

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::ImmersedFlowCaseConfig &config,
    const CheckpointV3ReadModules &modules, FlowState &state,
    const std::filesystem::path &directory) {
  return with_presence(
      read_checkpoint_v3_modules(mpi, decomposition, topology, geometry,
                                 boundaries, config, modules, state,
                                 directory),
      modules.presence);
}

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::ImmersedFlowCaseConfig &config,
    const immersed::ImmersedSurface &surface,
    const immersed::SurfaceQuery &query,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost,
    const immersed::WallQuadraturePlan &wall,
    const immersed::LocalFlowPatternTransform &transform,
    const FlowState &state, const FixedStepImmersedFlow &flow,
    CheckpointV3ControlState control,
    const std::filesystem::path &directory) {
  return write_checkpoint_v3(
      mpi, decomposition, topology, geometry, boundaries, config,
      {CheckpointV3Presence::constant_static_ibm, &surface, &query, &domain,
       &ghost, &wall, &transform, nullptr, &flow, nullptr},
      state, control, directory);
}

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::ImmersedFlowCaseConfig &config,
    const immersed::ImmersedSurface &surface,
    const immersed::SurfaceQuery &query,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost,
    const immersed::WallQuadraturePlan &wall,
    const immersed::LocalFlowPatternTransform &transform, FlowState &state,
    FixedStepImmersedFlow &flow, const std::filesystem::path &directory) {
  return read_checkpoint_v3(
      mpi, decomposition, topology, geometry, boundaries, config,
      {CheckpointV3Presence::constant_static_ibm, &surface, &query, &domain,
       &ghost, &wall, &transform, nullptr, &flow, nullptr},
      state, directory);
}

CheckpointV3Operation CheckpointV3Report::operation() const noexcept {
  return operation_;
}
CheckpointV3Disposition CheckpointV3Report::disposition() const noexcept {
  return disposition_;
}
CheckpointV3FailureReason CheckpointV3Report::reason() const noexcept {
  return reason_;
}
CheckpointV3Phase CheckpointV3Report::phase() const noexcept { return phase_; }
CheckpointV3Presence CheckpointV3Report::presence() const noexcept {
  return presence_;
}

diagnostics::Stage3PerformanceCounters
checkpoint_v3_performance_counters() noexcept {
  diagnostics::Stage3PerformanceCounters result;
  result.checkpoint_logical_io_bytes =
      checkpoint_logical_io_bytes.load(std::memory_order_relaxed);
  return result;
}
int CheckpointV3Report::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}
std::uint64_t CheckpointV3Report::step() const noexcept { return step_; }
double CheckpointV3Report::time_s() const noexcept { return time_s_; }
std::uint64_t CheckpointV3Report::manifest_crc64() const noexcept {
  return manifest_crc64_;
}
CheckpointV3CheckStatus CheckpointV3Report::crc_status() const noexcept {
  return crc_status_;
}
CheckpointV3CheckStatus
CheckpointV3Report::fingerprint_status() const noexcept {
  return fingerprint_status_;
}
CheckpointV3CheckStatus
CheckpointV3Report::partition_status() const noexcept {
  return partition_status_;
}
CheckpointV3CheckStatus CheckpointV3Report::rollback_status() const noexcept {
  return rollback_status_;
}

const CheckpointV3Report &CheckpointV3ReadResult::report() const noexcept {
  return report_;
}
bool CheckpointV3ReadResult::restored() const noexcept { return restored_; }
const CheckpointV3ControlState &CheckpointV3ReadResult::control_state() const {
  if (!restored_) {
    throw runtime::Error("Checkpoint v3 did not restore control state");
  }
  return control_;
}

} // namespace hundun::flow

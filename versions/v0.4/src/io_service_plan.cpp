// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_io.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kIoServicePlan = 10101U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

}  // namespace

Status IoServicePlan::compile(
    Span<const SnapshotFieldSpec> snapshot_fields,
    Span<const RuntimeServiceCapacity> services,
    std::size_t local_cells, IoServicePlan& out) noexcept {
  if (snapshot_fields.data == nullptr || snapshot_fields.size == 0U ||
      snapshot_fields.size > kMaximumSnapshotFields ||
      services.data == nullptr || services.size != kMaximumServices ||
      local_cells == 0U || out.fingerprint_ != 0U) {
    return {StatusCode::invalid_plan, kIoServicePlan};
  }
  std::size_t snapshot_doubles = 0U;
  for (std::size_t index = 0U; index < snapshot_fields.size; ++index) {
    const SnapshotFieldSpec field = snapshot_fields.data[index];
    if (field.components == 0U) {
      return {StatusCode::invalid_plan, kIoServicePlan};
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (snapshot_fields.data[prior].field == field.field) {
        return {StatusCode::invalid_plan, kIoServicePlan};
      }
    }
    if (field.components >
            std::numeric_limits<std::size_t>::max() - snapshot_doubles) {
      return {StatusCode::invalid_plan, kIoServicePlan};
    }
    snapshot_doubles += field.components;
  }
  if (snapshot_doubles > std::numeric_limits<std::size_t>::max() /
                             local_cells ||
      snapshot_doubles * local_cells >
          std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    return {StatusCode::invalid_plan, kIoServicePlan};
  }
  const std::size_t snapshot_bytes =
      snapshot_doubles * local_cells * sizeof(double);
  std::array<std::uint8_t, kMaximumServices> seen{};
  std::size_t maximum_staging = 0U;
  for (std::size_t index = 0U; index < services.size; ++index) {
    const RuntimeServiceCapacity service = services.data[index];
    const std::size_t kind = static_cast<std::size_t>(service.kind);
    if (kind >= seen.size() || seen[kind] != 0U || service.stage == 0U ||
        service.maximum_snapshot_bytes_per_rank < snapshot_bytes ||
        service.maximum_staging_bytes_per_rank == 0U ||
        service.maximum_collectives == 0U) {
      return {StatusCode::invalid_plan, kIoServicePlan};
    }
    seen[kind] = 1U;
    maximum_staging =
        std::max(maximum_staging, service.maximum_staging_bytes_per_rank);
  }
  if (std::any_of(seen.begin(), seen.end(),
                  [](std::uint8_t value) { return value == 0U; })) {
    return {StatusCode::invalid_plan, kIoServicePlan};
  }
  IoServicePlan candidate;
  std::copy(snapshot_fields.data,
            snapshot_fields.data + snapshot_fields.size,
            candidate.snapshot_fields_.begin());
  std::copy(services.data, services.data + services.size,
            candidate.services_.begin());
  candidate.snapshot_field_count_ = snapshot_fields.size;
  candidate.service_count_ = services.size;
  candidate.maximum_staging_bytes_ = maximum_staging;
  std::uint64_t hash = kFnvOffset;
  hash = mix(hash, local_cells);
  hash = mix(hash, snapshot_bytes);
  for (std::size_t index = 0U; index < snapshot_fields.size; ++index) {
    hash = mix(hash, snapshot_fields.data[index].field);
    hash = mix(hash, snapshot_fields.data[index].components);
  }
  for (std::size_t index = 0U; index < services.size; ++index) {
    const RuntimeServiceCapacity service = services.data[index];
    hash = mix(hash, static_cast<std::uint8_t>(service.kind));
    hash = mix(hash, service.stage);
    hash = mix(hash, service.maximum_snapshot_bytes_per_rank);
    hash = mix(hash, service.maximum_staging_bytes_per_rank);
    hash = mix(hash, service.maximum_collectives);
  }
  candidate.fingerprint_ = hash == 0U ? 1U : hash;
  out = candidate;
  return {};
}

}  // namespace hundun::v04

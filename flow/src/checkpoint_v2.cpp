// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/checkpoint_v2.hpp"
#include "adaptive_time_control_detail.hpp"
#include "checkpoint_v2_detail.hpp"
#include "checkpoint_v2_protocol.hpp"

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

template <class T>
constexpr std::uint8_t value(T item) noexcept {
  return static_cast<std::uint8_t>(item);
}

std::uint64_t report_seal(const CheckpointV2Report &report) noexcept {
  runtime::checkpoint_v2::Encoder encoder;
  const std::string domain = "hundun.checkpoint-v2.report-semantic.v1";
  encoder.string(domain);
  encoder.u8(value(report.operation()));
  encoder.u8(value(report.disposition()));
  encoder.u8(value(report.reason()));
  encoder.u8(value(report.phase()));
  encoder.i32(report.rank());
  encoder.i32(report.lowest_failing_rank());
  encoder.u64(report.step());
  encoder.f64(report.time_s());
  encoder.u64(report.local_logical_bytes());
  encoder.u64(report.local_actual_bytes());
  encoder.u64(report.global_logical_bytes());
  encoder.u64(report.global_actual_bytes());
  encoder.u64(report.local_crc64());
  encoder.u64(report.manifest_crc64());
  encoder.u64(report.file_count());
  encoder.u64(report.crc_check_count());
  encoder.u64(report.collective_count());
  encoder.u8(value(report.rank_crc_status()));
  encoder.u8(value(report.manifest_crc_status()));
  encoder.u8(value(report.exact_size_and_eof_status()));
  encoder.u8(value(report.fingerprint_status()));
  encoder.u8(value(report.partition_status()));
  encoder.u8(value(report.transaction_entry_status()));
  encoder.u8(value(report.publication_status()));
  encoder.u8(value(report.rollback_status()));
  return runtime::checkpoint_v2::crc64_ecma(encoder.bytes().data(),
                                            encoder.bytes().size());
}

using ByteVector = std::vector<std::uint8_t>;
using runtime::checkpoint_v2::Decoder;
using runtime::checkpoint_v2::Encoder;

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}
bool same(runtime::Box3 left, runtime::Box3 right) noexcept {
  return same(left.begin, right.begin) && same(left.end, right.end);
}
bool fp_equal(double left, double right) noexcept {
  std::uint64_t l{};
  std::uint64_t r{};
  std::memcpy(&l, &left, sizeof(l));
  std::memcpy(&r, &right, sizeof(r));
  return l == r;
}

void append_optional(Encoder &encoder, const std::optional<double> &item) {
  encoder.boolean(item.has_value());
  if (item)
    encoder.f64(*item);
}
void append_optional(Encoder &encoder,
                     const std::optional<runtime::Real3> &item) {
  encoder.boolean(item.has_value());
  if (item) {
    encoder.f64(item->x);
    encoder.f64(item->y);
    encoder.f64(item->z);
  }
}

template <class Function>
std::uint64_t make_fingerprint(std::string_view domain, Function &&function) {
  Encoder encoder;
  encoder.raw(domain.data(), domain.size());
  encoder.u8(0U);
  encoder.u32(1U);
  function(encoder);
  return runtime::checkpoint_v2::crc64_ecma(encoder.bytes().data(),
                                             encoder.bytes().size());
}

std::uint64_t resolved_fingerprint(const config::FlowCaseConfig &config) {
  return make_fingerprint("hundun.checkpoint-v2.resolved-case.v1",
                          [&](Encoder &encoder) {
    encoder.i32(config.schema_version);
    encoder.u8(value(config.simulation_type));
    encoder.u8(value(config.density_model));
    encoder.boolean(config.resources.expected_ranks.has_value());
    if (config.resources.expected_ranks)
      encoder.i32(*config.resources.expected_ranks);
    encoder.boolean(config.resources.process_grid.has_value());
    if (config.resources.process_grid) {
      encoder.i32(config.resources.process_grid->x);
      encoder.i32(config.resources.process_grid->y);
      encoder.i32(config.resources.process_grid->z);
    }
    encoder.i32(config.mesh.cells.x);
    encoder.i32(config.mesh.cells.y);
    encoder.i32(config.mesh.cells.z);
    for (double item : {config.mesh.origin_m.x, config.mesh.origin_m.y,
                        config.mesh.origin_m.z, config.mesh.length_m.x,
                        config.mesh.length_m.y, config.mesh.length_m.z})
      encoder.f64(item);
    encoder.u8(value(config.mesh.mapping));
    append_optional(encoder, config.mesh.warp_amplitude);
    encoder.u8(value(config.time.mode));
    encoder.i32(config.time.steps);
    for (double item : {
             config.time.initial_dt_s, config.time.min_dt_s,
             config.time.max_dt_s, config.time.cfl_target,
             config.time.diffusion_number_target, config.time.growth_factor,
             config.time.retry_factor})
      encoder.f64(item);
    encoder.i32(config.time.max_retries);
    encoder.f64(config.physics.rho_ref_kg_per_m3);
    encoder.f64(config.physics.dynamic_viscosity_pa_s);
    encoder.f64(config.physics.inlet_consistency_rtol);
    append_optional(encoder, config.physics.cp_J_per_kg_K);
    append_optional(encoder, config.physics.gas_constant_J_per_kg_K);
    append_optional(encoder, config.physics.thermodynamic_pressure_pa);
    encoder.u64(config.scalars.size());
    for (const auto &scalar : config.scalars) {
      encoder.string(scalar.name);
      encoder.f64(scalar.diffusivity_m2_per_s);
    }
    encoder.u64(config.boundaries.size());
    for (const auto &boundary : config.boundaries) {
      encoder.u8(value(boundary.patch));
      encoder.u8(value(boundary.type));
      append_optional(encoder, boundary.velocity_m_per_s);
      encoder.boolean(boundary.thermal_authority.has_value());
      if (boundary.thermal_authority)
        encoder.u8(value(*boundary.thermal_authority));
      append_optional(encoder, boundary.temperature_K);
      append_optional(encoder, boundary.enthalpy_J_per_kg);
      append_optional(encoder, boundary.density_kg_per_m3);
      encoder.boolean(boundary.scalar_values.has_value());
      if (boundary.scalar_values) {
        encoder.u64(boundary.scalar_values->size());
        for (const auto &item : *boundary.scalar_values) {
          encoder.string(item.name);
          encoder.f64(item.value);
        }
      }
      append_optional(encoder, boundary.pressure_perturbation_pa);
    }
  });
}

std::uint64_t topology_fingerprint(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, int rank_count) {
  return make_fingerprint("hundun.checkpoint-v2.topology-common.v1",
                          [&](Encoder &encoder) {
    const auto global = topology.global_extent();
    const auto grid = decomposition.process_grid();
    encoder.i32(rank_count);
    encoder.i32(grid.x);
    encoder.i32(grid.y);
    encoder.i32(grid.z);
    encoder.i32(global.x);
    encoder.i32(global.y);
    encoder.i32(global.z);
    for (const bool periodic : decomposition.periodic())
      encoder.boolean(periodic);
    encoder.u64(topology.global_cell_count());
    encoder.u64(topology.global_face_count());
    for (const auto &patch : topology.patches()) {
      encoder.u32(patch.stable_id());
      encoder.string(std::string(patch.name()));
      encoder.u8(value(patch.pairing_kind()));
      encoder.boolean(patch.paired_patch_id().has_value());
      if (patch.paired_patch_id())
        encoder.u32(*patch.paired_patch_id());
    }
  });
}

std::uint64_t geometry_fingerprint(const config::FlowCaseConfig &config) {
  return make_fingerprint("hundun.checkpoint-v2.geometry-common.v1",
                          [&](Encoder &encoder) {
    encoder.u8(value(config.mesh.mapping));
    encoder.i32(config.mesh.cells.x);
    encoder.i32(config.mesh.cells.y);
    encoder.i32(config.mesh.cells.z);
    for (double item : {config.mesh.origin_m.x, config.mesh.origin_m.y,
                        config.mesh.origin_m.z, config.mesh.length_m.x,
                        config.mesh.length_m.y, config.mesh.length_m.z})
      encoder.f64(item);
    append_optional(encoder, config.mesh.warp_amplitude);
  });
}

std::uint64_t boundary_fingerprint(const config::FlowCaseConfig &config) {
  return make_fingerprint("hundun.checkpoint-v2.boundary.v1",
                          [&](Encoder &encoder) {
    encoder.u64(config.scalars.size());
    for (const auto &item : config.scalars)
      encoder.string(item.name);
    for (const auto &item : config.boundaries) {
      encoder.u8(value(item.patch));
      encoder.u8(value(item.type));
      append_optional(encoder, item.velocity_m_per_s);
      append_optional(encoder, item.temperature_K);
      append_optional(encoder, item.enthalpy_J_per_kg);
      append_optional(encoder, item.density_kg_per_m3);
      append_optional(encoder, item.pressure_perturbation_pa);
    }
  });
}

std::vector<runtime::FieldId> ordered_fields(const FlowFieldIds &fields) {
  std::vector<runtime::FieldId> result{
      fields.density, fields.velocity, fields.mechanical_pressure,
      fields.face_velocity, fields.face_mass_flux};
  result.insert(result.end(), fields.transported_cell_fields.begin(),
                fields.transported_cell_fields.end());
  return result;
}

std::uint64_t field_schema_fingerprint(const runtime::FieldRegistry &registry,
                                       const FlowFieldIds &fields) {
  return make_fingerprint("hundun.checkpoint-v2.field-schema.v1",
                          [&](Encoder &encoder) {
    const auto ids = ordered_fields(fields);
    encoder.u64(ids.size());
    for (std::size_t index = 0; index < ids.size(); ++index) {
      const auto id = ids[index];
      const auto &descriptor = registry.descriptor(id);
      const std::uint8_t role =
          index < 5U ? static_cast<std::uint8_t>(index) : 5U;
      encoder.u8(role);
      encoder.u32(index < 5U ? 0U
                             : static_cast<std::uint32_t>(index - 5U));
      encoder.u32(id);
      encoder.string(descriptor.name);
      encoder.string(descriptor.unit);
      encoder.string(descriptor.owner);
      encoder.u8(value(descriptor.space));
      encoder.u8(value(descriptor.scalar_type));
      encoder.u32(descriptor.components);
      encoder.i32(descriptor.ghost_width);
      encoder.boolean(descriptor.conservative);
      encoder.u8(value(descriptor.restart));
      encoder.u8(value(descriptor.output));
    }
  });
}

std::uint64_t local_layout_fingerprint(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, int rank) {
  return make_fingerprint("hundun.checkpoint-v2.local-layout.v1",
                          [&](Encoder &encoder) {
    const auto box = decomposition.owned_box();
    const auto extent = decomposition.local_extent();
    const auto grid = decomposition.process_grid();
    encoder.i32(rank);
    encoder.i32(grid.x);
    encoder.i32(grid.y);
    encoder.i32(grid.z);
    for (int item : {box.begin.x, box.begin.y, box.begin.z, box.end.x,
                     box.end.y, box.end.z, extent.x, extent.y, extent.z})
      encoder.i32(item);
    encoder.u64(topology.local_cell_count());
    for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
      encoder.u8(value(topology.cell_ownership(cell)));
      encoder.u64(topology.global_cell_id(cell));
    }
    encoder.u64(topology.local_face_count());
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      encoder.u8(value(topology.face_ownership(face)));
      encoder.u64(topology.global_face_id(face));
    }
  });
}

std::array<std::uint64_t, 5> fingerprints(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const config::FlowCaseConfig &config,
    const runtime::FieldRegistry &registry, const FlowFieldIds &fields,
    int ranks) {
  return {resolved_fingerprint(config),
          topology_fingerprint(decomposition, topology, ranks),
          geometry_fingerprint(config), boundary_fingerprint(config),
          field_schema_fingerprint(registry, fields)};
}

ByteVector encode_global_payload(
    AcceptedStepMetadata metadata, const TimeControlState &time,
    const std::optional<IdealGasClosureState> &closure) {
  Encoder encoder;
  encoder.u32(1U);
  encoder.u64(metadata.step);
  encoder.f64(metadata.time_s);
  encoder.f64(metadata.dt_s);
  encoder.f64(metadata.previous_dt_s);
  encoder.u8(value(metadata.order));
  encoder.u32(time.schema_version);
  encoder.u64(time.accepted_step);
  encoder.f64(time.proposed_next_dt_s);
  encoder.f64(time.last_accepted_dt_s);
  encoder.u8(value(time.last_accepted_order));
  encoder.boolean(time.history_ready);
  encoder.boolean(time.last_all_linear_solves_within_half_limit);
  encoder.f64(time.last_convective_rate_per_s);
  encoder.f64(time.last_diffusive_rate_per_s);
  encoder.boolean(time.last_stability_metrics_available);
  encoder.u32(time.last_retry_count);
  encoder.u64(time.revision);
  encoder.u64(time.state_seal);
  encoder.boolean(closure.has_value());
  if (closure) {
    encoder.u8(value(closure->mode));
    encoder.f64(closure->thermodynamic_pressure_pa);
    encoder.boolean(closure->target_mass_kg.has_value());
    if (closure->target_mass_kg)
      encoder.f64(*closure->target_mass_kg);
    encoder.u64(closure->revision);
  }
  return std::move(encoder).take();
}

struct GlobalPayload final {
  AcceptedStepMetadata metadata;
  TimeControlState time;
  std::optional<IdealGasClosureState> closure;
};

GlobalPayload decode_global_payload(const ByteVector &bytes) {
  Decoder decoder(bytes);
  if (decoder.u32() != 1U)
    throw runtime::Error("Checkpoint v2 global payload schema is invalid");
  GlobalPayload result;
  result.metadata.step = decoder.u64();
  result.metadata.time_s = decoder.f64();
  result.metadata.dt_s = decoder.f64();
  result.metadata.previous_dt_s = decoder.f64();
  const auto order = decoder.u8();
  if (order > 1U)
    throw runtime::Error("Checkpoint v2 time order is invalid");
  result.metadata.order = static_cast<MomentumTimeOrder>(order);
  result.time.schema_version = decoder.u32();
  result.time.accepted_step = decoder.u64();
  result.time.proposed_next_dt_s = decoder.f64();
  result.time.last_accepted_dt_s = decoder.f64();
  const auto last_order = decoder.u8();
  if (last_order > 1U)
    throw runtime::Error("Checkpoint v2 controller order is invalid");
  result.time.last_accepted_order =
      static_cast<MomentumTimeOrder>(last_order);
  result.time.history_ready = decoder.boolean();
  result.time.last_all_linear_solves_within_half_limit = decoder.boolean();
  result.time.last_convective_rate_per_s = decoder.f64();
  result.time.last_diffusive_rate_per_s = decoder.f64();
  result.time.last_stability_metrics_available = decoder.boolean();
  result.time.last_retry_count = decoder.u32();
  result.time.revision = decoder.u64();
  result.time.state_seal = decoder.u64();
  if (decoder.boolean()) {
    IdealGasClosureState closure;
    const auto mode = decoder.u8();
    if (mode > 1U)
      throw runtime::Error("Checkpoint v2 pressure mode is invalid");
    closure.mode = static_cast<IdealGasPressureMode>(mode);
    closure.thermodynamic_pressure_pa = decoder.f64();
    if (decoder.boolean())
      closure.target_mass_kg = decoder.f64();
    closure.revision = decoder.u64();
    result.closure = closure;
  }
  decoder.require_eof();
  return result;
}

void append_record(Encoder &encoder, std::uint8_t layer, std::uint8_t role,
                   std::uint32_t transported, runtime::FieldId field,
                   runtime::FunctionSpace space, std::uint32_t components,
                   const std::vector<double> &values,
                   std::uint64_t &logical_bytes) {
  if (values.size() >
      std::numeric_limits<std::uint64_t>::max() / sizeof(double))
    throw runtime::Error("Checkpoint v2 field value size overflows");
  const auto bytes =
      static_cast<std::uint64_t>(values.size() * sizeof(double));
  encoder.u8(layer);
  encoder.u8(role);
  encoder.u32(transported);
  encoder.u32(field);
  encoder.u8(space == runtime::FunctionSpace::cell_average ? 0U : 1U);
  encoder.u8(0U);
  encoder.u32(components);
  encoder.u64(values.size() / components);
  encoder.u64(bytes);
  for (double item : values)
    encoder.f64(item);
  if (logical_bytes > std::numeric_limits<std::uint64_t>::max() - bytes)
    throw runtime::Error("Checkpoint v2 logical byte count overflows");
  logical_bytes += bytes;
}

ByteVector encode_rank_payload(const FlowState &state,
                               std::uint64_t &logical_bytes) {
  const auto extent =
      detail::FlowStateCheckpointAccess::layout(state).cell_interior_extent;
  const auto face_count =
      detail::FlowStateCheckpointAccess::layout(state).face_count;
  const auto history = state.snapshot(FlowLayer::history);
  const auto committed = state.snapshot(FlowLayer::committed);
  const auto &fields = state.fields();
  Encoder encoder;
  encoder.u32(1U);
  encoder.i32(extent.x);
  encoder.i32(extent.y);
  encoder.i32(extent.z);
  encoder.u64(face_count);
  encoder.u32(static_cast<std::uint32_t>(
      fields.transported_cell_fields.size()));
  encoder.u32(static_cast<std::uint32_t>(
      2U * (5U + fields.transported_cell_fields.size())));
  logical_bytes = 0U;
  const auto layer = [&](std::uint8_t layer_id,
                         const FlowLayerValues &values) {
    append_record(encoder, layer_id, 0U, 0U, fields.density,
                  runtime::FunctionSpace::cell_average, 1U, values.density,
                  logical_bytes);
    append_record(encoder, layer_id, 1U, 0U, fields.velocity,
                  runtime::FunctionSpace::cell_average, 3U, values.velocity,
                  logical_bytes);
    append_record(encoder, layer_id, 2U, 0U, fields.mechanical_pressure,
                  runtime::FunctionSpace::cell_average, 1U,
                  values.mechanical_pressure, logical_bytes);
    append_record(encoder, layer_id, 3U, 0U, fields.face_velocity,
                  runtime::FunctionSpace::face_value, 3U,
                  values.face_velocity, logical_bytes);
    append_record(encoder, layer_id, 4U, 0U, fields.face_mass_flux,
                  runtime::FunctionSpace::face_value, 1U,
                  values.face_mass_flux, logical_bytes);
    for (std::size_t index = 0;
         index < fields.transported_cell_fields.size(); ++index)
      append_record(
          encoder, layer_id, 5U, static_cast<std::uint32_t>(index),
          fields.transported_cell_fields[index],
          runtime::FunctionSpace::cell_average, 1U,
          values.transported_cell_fields[index], logical_bytes);
  };
  layer(0U, history);
  layer(1U, committed);
  return std::move(encoder).take();
}

struct DecodedRecord final {
  std::uint8_t layer{};
  std::uint8_t role{};
  std::uint32_t transported{};
  runtime::FieldId field{};
  std::uint8_t space{};
  std::uint32_t components{};
  std::vector<double> values;
};

std::pair<FlowLayerValues, FlowLayerValues>
decode_rank_payload(const ByteVector &bytes, const FlowState &state) {
  Decoder decoder(bytes);
  const auto layout = detail::FlowStateCheckpointAccess::layout(state);
  const auto &fields = state.fields();
  if (decoder.u32() != 1U ||
      !same(runtime::Int3{decoder.i32(), decoder.i32(), decoder.i32()},
            layout.cell_interior_extent) ||
      decoder.u64() != layout.face_count ||
      decoder.u32() != fields.transported_cell_fields.size() ||
      decoder.u32() != 2U * (5U + fields.transported_cell_fields.size()))
    throw runtime::Error("Checkpoint v2 rank payload layout is invalid");
  FlowLayerValues layers[2];
  const std::size_t cell_count =
      static_cast<std::size_t>(layout.cell_interior_extent.x) *
      static_cast<std::size_t>(layout.cell_interior_extent.y) *
      static_cast<std::size_t>(layout.cell_interior_extent.z);
  for (std::size_t layer = 0; layer < 2U; ++layer) {
    for (std::size_t role = 0;
         role < 5U + fields.transported_cell_fields.size(); ++role) {
      const auto observed_layer = decoder.u8();
      const auto observed_role = decoder.u8();
      const auto transported = decoder.u32();
      const auto field = decoder.u32();
      const auto space = decoder.u8();
      const auto scalar = decoder.u8();
      const auto components = decoder.u32();
      const auto items = decoder.u64();
      const auto byte_size = decoder.u64();
      const std::uint8_t expected_role =
          role < 5U ? static_cast<std::uint8_t>(role) : 5U;
      const std::uint32_t expected_transport =
          role < 5U ? 0U : static_cast<std::uint32_t>(role - 5U);
      const auto expected_field =
          role < 5U ? ordered_fields(fields)[role]
                    : fields.transported_cell_fields[role - 5U];
      const std::uint32_t expected_components =
          role == 1U || role == 3U ? 3U : 1U;
      const bool face = role == 3U || role == 4U;
      const std::size_t expected_items = face ? layout.face_count : cell_count;
      if (observed_layer != layer || observed_role != expected_role ||
          transported != expected_transport || field != expected_field ||
          space != (face ? 1U : 0U) || scalar != 0U ||
          components != expected_components || items != expected_items ||
          byte_size != items * components * sizeof(double))
        throw runtime::Error("Checkpoint v2 field record is invalid");
      std::vector<double> values;
      values.reserve(static_cast<std::size_t>(items) * components);
      for (std::uint64_t index = 0; index < items * components; ++index)
        values.push_back(decoder.f64());
      auto &target = layers[layer];
      switch (role) {
      case 0:
        target.density = std::move(values);
        break;
      case 1:
        target.velocity = std::move(values);
        break;
      case 2:
        target.mechanical_pressure = std::move(values);
        break;
      case 3:
        target.face_velocity = std::move(values);
        break;
      case 4:
        target.face_mass_flux = std::move(values);
        break;
      default:
        target.transported_cell_fields.push_back(std::move(values));
        break;
      }
    }
  }
  decoder.require_eof();
  return {std::move(layers[0]), std::move(layers[1])};
}

bool valid_layer(const FlowLayerValues &layer,
                 const config::FlowCaseConfig &config) {
  const auto finite = [](const auto &values) {
    return std::all_of(values.begin(), values.end(),
                       [](double item) { return std::isfinite(item); });
  };
  if (!finite(layer.density) || !finite(layer.velocity) ||
      !finite(layer.mechanical_pressure) || !finite(layer.face_velocity) ||
      !finite(layer.face_mass_flux) ||
      !std::all_of(layer.density.begin(), layer.density.end(),
                   [](double item) { return item > 0.0; }))
    return false;
  for (const auto &transported : layer.transported_cell_fields)
    if (!finite(transported))
      return false;
  if (config.density_model == config::DensityModel::constant)
    for (double item : layer.density)
      if (!fp_equal(item, config.physics.rho_ref_kg_per_m3))
        return false;
  if (config.density_model == config::DensityModel::ideal_gas) {
    if (layer.transported_cell_fields.empty() ||
        !std::all_of(layer.transported_cell_fields.front().begin(),
                     layer.transported_cell_fields.front().end(),
                     [](double item) {
                       return item > 0.0 && std::isfinite(item);
                     }))
      return false;
  }
  return true;
}

double relative_error(double observed, double reference) noexcept {
  return std::abs(observed - reference) /
         std::max(std::abs(reference), std::numeric_limits<double>::min());
}

bool ideal_gas_layers_valid(
    const runtime::MpiContext &mpi, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::FlowCaseConfig &config, const FlowLayerValues &history,
    const FlowLayerValues &committed, const IdealGasClosureState &closure,
    std::uint64_t &collective_count) {
  if (!config.physics.cp_J_per_kg_K ||
      !config.physics.gas_constant_J_per_kg_K ||
      !config.physics.thermodynamic_pressure_pa ||
      !(*config.physics.cp_J_per_kg_K > 0.0) ||
      !(*config.physics.gas_constant_J_per_kg_K > 0.0) ||
      !(*config.physics.thermodynamic_pressure_pa > 0.0))
    return false;
  const double cp = *config.physics.cp_J_per_kg_K;
  const double gas_constant = *config.physics.gas_constant_J_per_kg_K;
  const bool open = boundaries.open_domain();
  bool valid =
      closure.mode == (open ? IdealGasPressureMode::open_fixed
                            : IdealGasPressureMode::closed_dynamic) &&
      closure.thermodynamic_pressure_pa > 0.0 &&
      std::isfinite(closure.thermodynamic_pressure_pa) &&
      closure.revision != std::numeric_limits<std::uint64_t>::max() &&
      closure.target_mass_kg.has_value() == !open;
  if (open)
    valid =
        valid &&
        fp_equal(closure.thermodynamic_pressure_pa,
                 *config.physics.thermodynamic_pressure_pa);
  const double target_mass =
      closure.target_mass_kg.value_or(
          std::numeric_limits<double>::quiet_NaN());
  if (!open)
    valid = valid && target_mass > 0.0 && std::isfinite(target_mass);

  std::array<double, 3> local_sums{};
  bool layer_shapes_valid = true;
  const auto accumulate = [&](const FlowLayerValues &layer,
                              std::size_t mass_index,
                              bool accumulate_inverse_temperature) {
    if (layer.transported_cell_fields.empty() ||
        layer.density.size() != topology.owned_cell_count() ||
        layer.transported_cell_fields.front().size() !=
            topology.owned_cell_count()) {
      valid = false;
      layer_shapes_valid = false;
      return;
    }
    for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double rho = layer.density[cell];
      const double q = layer.transported_cell_fields.front()[cell];
      const double enthalpy = q / rho;
      const double temperature = enthalpy / cp;
      const double volume = geometry.cell_volume_m3(cell);
      if (!(rho > 0.0) || !(q > 0.0) || !(enthalpy > 0.0) ||
          !(temperature > 0.0) || !(volume > 0.0) ||
          !std::isfinite(rho) || !std::isfinite(q) ||
          !std::isfinite(enthalpy) || !std::isfinite(temperature) ||
          !std::isfinite(volume)) {
        valid = false;
        continue;
      }
      local_sums[mass_index] += volume * rho;
      if (accumulate_inverse_temperature)
        local_sums[2] += volume / temperature;
    }
  };
  accumulate(history, 0U, !open);
  accumulate(committed, 1U, false);
  if (!open) {
    mpi.allreduce_fp64_in_place(local_sums.data(), local_sums.size(),
                                runtime::Fp64ReductionOperation::sum);
    ++collective_count;
    valid = valid && local_sums[2] > 0.0 &&
            std::isfinite(local_sums[2]) &&
            relative_error(local_sums[0], target_mass) <= 5.0e-12 &&
            relative_error(local_sums[1], target_mass) <= 5.0e-12;
  }
  if (!layer_shapes_valid)
    return false;

  const double history_pressure =
      open ? closure.thermodynamic_pressure_pa
           : target_mass * gas_constant / local_sums[2];
  const auto eos_valid = [&](const FlowLayerValues &layer, double pressure) {
    if (!(pressure > 0.0) || !std::isfinite(pressure))
      return false;
    bool result = true;
    for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double rho = layer.density[cell];
      const double q = layer.transported_cell_fields.front()[cell];
      const double temperature = (q / rho) / cp;
      const double ratio = rho * gas_constant * temperature / pressure;
      result =
          result && std::isfinite(ratio) &&
          std::abs(ratio - 1.0) <= 1.0e-12;
    }
    return result;
  };
  valid = valid && eos_valid(history, history_pressure) &&
          eos_valid(committed, closure.thermodynamic_pressure_pa);
  return valid;
}

std::string rank_filename(int rank) {
  std::array<char, 32> buffer{};
  const int count = std::snprintf(buffer.data(), buffer.size(),
                                  "rank-%06d.v2.bin", rank);
  if (count <= 0 || static_cast<std::size_t>(count) >= buffer.size())
    throw runtime::Error("Checkpoint v2 rank filename is invalid");
  return std::string(buffer.data(), static_cast<std::size_t>(count));
}

struct RankRecordWire final {
  std::int32_t rank{};
  std::int32_t begin[3]{};
  std::int32_t end[3]{};
  std::uint64_t logical{};
  std::uint64_t actual{};
  std::uint64_t crc{};
  std::uint64_t layout{};
};

bool path_agrees(const runtime::MpiContext &mpi,
                 const std::filesystem::path &path,
                 std::uint64_t &collective_count) {
  const auto text = path.lexically_normal().generic_string();
  const auto hash = make_fingerprint("hundun.checkpoint-v2.path.v1",
                                     [&](Encoder &encoder) {
    encoder.string(text);
  });
  std::array<std::uint64_t, 2> local{
      static_cast<std::uint64_t>(text.size()), hash};
  std::array<std::uint64_t, 2> minima{};
  std::array<std::uint64_t, 2> maxima{};
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), minima.data(), 2, MPI_UINT64_T, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(Checkpoint v2 path minimum)");
  ++collective_count;
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), maxima.data(), 2, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce(Checkpoint v2 path maximum)");
  ++collective_count;
  return !text.empty() && text.find('\0') == std::string::npos &&
         minima == maxima;
}

struct Convergence final {
  bool ok{};
  int failing_rank{-1};
};

Convergence converge(const runtime::MpiContext &mpi, bool local_ok,
                     std::uint64_t &collective_count,
                     std::string_view operation) {
  const int candidate = local_ok ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  ++collective_count;
  return {lowest == mpi.size(), lowest == mpi.size() ? -1 : lowest};
}

std::uint64_t common_fingerprint(
    const std::array<std::uint64_t, 5> &items, int ranks, runtime::Int3 grid,
    const ByteVector &global) {
  return make_fingerprint("hundun.checkpoint-v2.common.v1",
                          [&](Encoder &encoder) {
    for (auto item : items)
      encoder.u64(item);
    encoder.i32(ranks);
    encoder.i32(grid.x);
    encoder.i32(grid.y);
    encoder.i32(grid.z);
    encoder.u64(global.size());
    encoder.raw(global.data(), global.size());
  });
}

bool directory_inventory_valid(const std::filesystem::path &directory,
                               int ranks) {
  std::error_code error;
  const auto status = std::filesystem::symlink_status(directory, error);
  if (error || status.type() != std::filesystem::file_type::directory)
    return false;
  std::set<std::string> expected{"COMPLETED", "manifest.v2.bin"};
  for (int rank = 0; rank < ranks; ++rank)
    expected.insert(rank_filename(rank));
  std::set<std::string> observed;
  for (std::filesystem::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const auto entry_status = iterator->symlink_status(error);
    if (error || entry_status.type() != std::filesystem::file_type::regular)
      return false;
    observed.insert(iterator->path().filename().generic_string());
  }
  return !error && observed == expected;
}

} // namespace

#define HUNDUN_CHECKPOINT_GETTER(type, name)                                  \
  type CheckpointV2Report::name() const noexcept { return name##_; }
HUNDUN_CHECKPOINT_GETTER(CheckpointV2Operation, operation)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2Disposition, disposition)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2FailureReason, reason)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2Phase, phase)
HUNDUN_CHECKPOINT_GETTER(int, rank)
HUNDUN_CHECKPOINT_GETTER(int, lowest_failing_rank)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, step)
HUNDUN_CHECKPOINT_GETTER(double, time_s)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, local_logical_bytes)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, local_actual_bytes)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, global_logical_bytes)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, global_actual_bytes)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, local_crc64)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, manifest_crc64)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, file_count)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, crc_check_count)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, collective_count)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, rank_crc_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, manifest_crc_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, exact_size_and_eof_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, fingerprint_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, partition_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, transaction_entry_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, publication_status)
HUNDUN_CHECKPOINT_GETTER(CheckpointV2CheckStatus, rollback_status)
HUNDUN_CHECKPOINT_GETTER(std::uint64_t, semantic_fingerprint)
#undef HUNDUN_CHECKPOINT_GETTER

const CheckpointV2Report &CheckpointV2ReadResult::report() const noexcept {
  return report_;
}
CheckpointV2ReadResult::CheckpointV2ReadResult(CheckpointV2Report report)
    : report_(std::move(report)) {}
bool CheckpointV2ReadResult::restored() const noexcept { return restored_; }
const TimeControlState &CheckpointV2ReadResult::time_control_state() const {
  if (!restored_)
    throw runtime::Error("Checkpoint v2 did not restore time-control state");
  return time_control_;
}
bool CheckpointV2ReadResult::ideal_gas_closure_state_available() const noexcept {
  return restored_ && closure_.has_value();
}
const IdealGasClosureState &
CheckpointV2ReadResult::ideal_gas_closure_state() const {
  if (!ideal_gas_closure_state_available())
    throw runtime::Error("Checkpoint v2 has no ideal-gas closure state");
  return *closure_;
}

CheckpointV2Report detail::CheckpointV2Access::failed(
    CheckpointV2Operation operation, int rank, CheckpointV2FailureReason reason,
    CheckpointV2Phase phase) {
  CheckpointV2Report report;
  report.operation_ = operation;
  report.rank_ = rank;
  report.reason_ = reason;
  report.phase_ = phase;
  report.lowest_failing_rank_ = rank;
  authenticate(report);
  return report;
}
CheckpointV2Report
detail::CheckpointV2Access::make(CheckpointV2ReportValues values) {
  CheckpointV2Report report;
  report.operation_ = values.operation;
  report.disposition_ = values.disposition;
  report.reason_ = values.reason;
  report.phase_ = values.phase;
  report.rank_ = values.rank;
  report.lowest_failing_rank_ = values.lowest_failing_rank;
  report.step_ = values.step;
  report.time_s_ = values.time_s;
  report.local_logical_bytes_ = values.local_logical_bytes;
  report.local_actual_bytes_ = values.local_actual_bytes;
  report.global_logical_bytes_ = values.global_logical_bytes;
  report.global_actual_bytes_ = values.global_actual_bytes;
  report.local_crc64_ = values.local_crc64;
  report.manifest_crc64_ = values.manifest_crc64;
  report.file_count_ = values.file_count;
  report.crc_check_count_ = values.crc_check_count;
  report.collective_count_ = values.collective_count;
  report.rank_crc_status_ = values.rank_crc;
  report.manifest_crc_status_ = values.manifest_crc;
  report.exact_size_and_eof_status_ = values.exact_size_eof;
  report.fingerprint_status_ = values.fingerprint;
  report.partition_status_ = values.partition;
  report.transaction_entry_status_ = values.transaction_entry;
  report.publication_status_ = values.publication;
  report.rollback_status_ = values.rollback;
  authenticate(report);
  return report;
}
CheckpointV2ReadResult detail::CheckpointV2Access::make_read(
    CheckpointV2Report report, TimeControlState time,
    std::optional<IdealGasClosureState> closure, bool restored) {
  CheckpointV2ReadResult result(std::move(report));
  result.restored_ = restored;
  result.time_control_ = time;
  result.closure_ = std::move(closure);
  return result;
}
void detail::CheckpointV2Access::authenticate(
    CheckpointV2Report &report) noexcept {
  report.semantic_fingerprint_ = report_seal(report);
}
CheckpointV2ReadResult detail::CheckpointV2Access::failed_read(
    int rank, CheckpointV2FailureReason reason, CheckpointV2Phase phase) {
  return CheckpointV2ReadResult(
      failed(CheckpointV2Operation::read, rank, reason, phase));
}

CheckpointV2Report write_checkpoint_v2(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::FlowCaseConfig &config, const FlowState &state,
    const TimeControlState &time,
    std::optional<IdealGasClosureState> closure,
    const std::filesystem::path &directory) {
  static_cast<void>(boundaries);
  detail::CheckpointV2ReportValues values;
  values.operation = CheckpointV2Operation::write;
  values.rank = mpi.rank();
  std::uint64_t collectives = 0U;
  const auto fail = [&](CheckpointV2FailureReason reason,
                        CheckpointV2Phase phase, int rank) {
    values.disposition = CheckpointV2Disposition::failed;
    values.reason = reason;
    values.phase = phase;
    values.lowest_failing_rank = rank;
    values.collective_count = collectives;
    return detail::CheckpointV2Access::make(values);
  };

  AcceptedStepMetadata metadata{};
  ByteVector global_payload;
  ByteVector rank_payload;
  ByteVector rank_wrapper;
  FlowLayerValues history;
  FlowLayerValues committed;
  std::array<std::uint64_t, 5> identity{};
  std::uint64_t local_layout{};
  runtime::Box3 owned{};
  bool local_ok = true;
  CheckpointV2FailureReason local_reason{
      CheckpointV2FailureReason::invalid_input};
  bool common_path = false;
  try {
    common_path = path_agrees(mpi, directory, collectives);
  } catch (...) {
    throw;
  }
  try {
    local_ok = detail::FlowStateCheckpointAccess::live(state) &&
               !detail::FlowStateCheckpointAccess::attempt_active(state) &&
               geometry.compatible(topology) &&
               same(topology.global_extent(), decomposition.global_extent()) &&
               same(topology.owned_global_box(), decomposition.owned_box()) &&
               topology.local_face_count() ==
                   detail::FlowStateCheckpointAccess::layout(state).face_count &&
               same(decomposition.local_extent(),
                    detail::FlowStateCheckpointAccess::layout(state)
                        .cell_interior_extent) &&
               config.schema_version == 2 &&
               config.simulation_type ==
                   config::SimulationType::variable_density_flow &&
               same(config.mesh.cells, topology.global_extent());
    local_ok = local_ok && common_path;
    if (!local_ok)
      throw runtime::Error("Checkpoint v2 preflight identity is invalid");
    metadata = state.metadata();
    if (!detail::TimeControlStateCodec::semantically_valid(
            config.time, config.density_model, metadata, time))
      throw runtime::Error("Checkpoint v2 time-control state is invalid");
    const bool ideal =
        config.density_model == config::DensityModel::ideal_gas;
    if (closure.has_value() != ideal)
      throw runtime::Error("Checkpoint v2 closure presence is invalid");
    if (closure &&
        (!(closure->thermodynamic_pressure_pa > 0.0) ||
         !std::isfinite(closure->thermodynamic_pressure_pa) ||
         closure->revision == std::numeric_limits<std::uint64_t>::max() ||
         (closure->mode == IdealGasPressureMode::closed_dynamic) !=
             closure->target_mass_kg.has_value() ||
         (closure->target_mass_kg &&
          (!(*closure->target_mass_kg > 0.0) ||
           !std::isfinite(*closure->target_mass_kg)))))
      throw runtime::Error("Checkpoint v2 closure state is invalid");
    history = state.snapshot(FlowLayer::history);
    committed = state.snapshot(FlowLayer::committed);
    if (!valid_layer(history, config) || !valid_layer(committed, config)) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 physical state is invalid");
    }
    const auto &registry = detail::FlowStateCheckpointAccess::registry(state);
    identity = fingerprints(decomposition, topology, config, registry,
                            state.fields(), mpi.size());
    local_layout =
        local_layout_fingerprint(decomposition, topology, mpi.rank());
    global_payload = encode_global_payload(metadata, time, closure);
    rank_payload = encode_rank_payload(state, values.local_logical_bytes);
    rank_wrapper = runtime::checkpoint_v2::encode_rank_wrapper(
        mpi.rank(), mpi.size(), rank_payload);
    values.local_actual_bytes = rank_wrapper.size();
    values.local_crc64 = runtime::checkpoint_v2::crc64_ecma(
        rank_wrapper.data(), rank_wrapper.size());
    owned = decomposition.owned_box();
  } catch (...) {
    local_ok = false;
  }
  auto status =
      converge(mpi, local_ok, collectives, "MPI_Allreduce(Checkpoint preflight)");
  if (!status.ok) {
    values.step = local_ok ? metadata.step : 0U;
    values.time_s = local_ok ? metadata.time_s : 0.0;
    values.fingerprint = CheckpointV2CheckStatus::not_checked;
    values.partition = local_reason == CheckpointV2FailureReason::layout
                           ? CheckpointV2CheckStatus::failed
                           : CheckpointV2CheckStatus::not_checked;
    return fail(local_reason, CheckpointV2Phase::preflight,
                status.failing_rank);
  }
  values.step = metadata.step;
  values.time_s = metadata.time_s;
  values.fingerprint = CheckpointV2CheckStatus::passed;
  values.partition = CheckpointV2CheckStatus::passed;

  if (config.density_model == config::DensityModel::ideal_gas) {
    const bool model_valid = ideal_gas_layers_valid(
        mpi, topology, geometry, boundaries, config, history, committed,
        *closure, collectives);
    status = converge(mpi, model_valid, collectives,
                      "MPI_Allreduce(Checkpoint ideal-gas write state)");
    if (!status.ok)
      return fail(CheckpointV2FailureReason::state,
                  CheckpointV2Phase::rank_payload, status.failing_rank);
  }

  bool filesystem_ok = true;
  if (mpi.rank() == 0) {
    try {
      std::error_code error;
      const auto target_status = std::filesystem::symlink_status(directory, error);
      const bool absent =
          error == std::errc::no_such_file_or_directory ||
          target_status.type() == std::filesystem::file_type::not_found;
      if (!absent || std::filesystem::is_symlink(target_status))
        throw runtime::Error("Checkpoint v2 target already exists");
      error.clear();
      auto parent = directory.parent_path();
      if (parent.empty())
        parent = ".";
      const auto parent_status = std::filesystem::symlink_status(parent, error);
      if (error ||
          parent_status.type() != std::filesystem::file_type::directory)
        throw runtime::Error("Checkpoint v2 parent is not a directory");
      if (!std::filesystem::create_directory(directory, error) || error)
        throw runtime::Error("Checkpoint v2 directory creation failed");
    } catch (...) {
      filesystem_ok = false;
    }
  }
  status = converge(mpi, filesystem_ok, collectives,
                    "MPI_Allreduce(Checkpoint directory)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::rank_temporary_file, status.failing_rank);

  const auto rank_name = rank_filename(mpi.rank());
  const auto rank_temp = directory / (rank_name + ".tmp");
  bool rank_written = true;
  try {
    const auto verified =
        runtime::checkpoint_v2::write_verified_temporary(rank_temp,
                                                          rank_wrapper);
    rank_written = verified.actual_size == rank_wrapper.size() &&
                   verified.crc64 == values.local_crc64;
  } catch (...) {
    rank_written = false;
  }
  status = converge(mpi, rank_written, collectives,
                    "MPI_Allreduce(Checkpoint rank temporary)");
  if (!status.ok) {
    values.rank_crc = rank_written ? CheckpointV2CheckStatus::passed
                                   : CheckpointV2CheckStatus::failed;
    values.exact_size_eof = rank_written ? CheckpointV2CheckStatus::passed
                                         : CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::rank_temporary_file, status.failing_rank);
  }
  values.rank_crc = CheckpointV2CheckStatus::passed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  ++values.crc_check_count;

  bool rank_published = true;
  try {
    runtime::checkpoint_v2::publish_no_overwrite(rank_temp,
                                                  directory / rank_name);
  } catch (...) {
    rank_published = false;
  }
  status = converge(mpi, rank_published, collectives,
                    "MPI_Allreduce(Checkpoint rank publish)");
  if (!status.ok) {
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::rank_publish, status.failing_rank);
  }

  RankRecordWire local_record;
  local_record.rank = mpi.rank();
  local_record.begin[0] = owned.begin.x;
  local_record.begin[1] = owned.begin.y;
  local_record.begin[2] = owned.begin.z;
  local_record.end[0] = owned.end.x;
  local_record.end[1] = owned.end.y;
  local_record.end[2] = owned.end.z;
  local_record.logical = values.local_logical_bytes;
  local_record.actual = values.local_actual_bytes;
  local_record.crc = values.local_crc64;
  local_record.layout = local_layout;
  std::vector<RankRecordWire> records;
  if (mpi.rank() == 0)
    records.resize(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(
      MPI_Gather(&local_record, static_cast<int>(sizeof(local_record)),
                 MPI_BYTE, records.data(), static_cast<int>(sizeof(local_record)),
                 MPI_BYTE, 0, mpi.comm()),
      "MPI_Gather(Checkpoint rank records)");
  ++collectives;

  std::uint64_t global_logical = 0U;
  std::uint64_t rank_actual = 0U;
  runtime::check_mpi_result(
      MPI_Allreduce(&values.local_logical_bytes, &global_logical, 1,
                    MPI_UINT64_T, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Checkpoint logical bytes)");
  ++collectives;
  runtime::check_mpi_result(
      MPI_Allreduce(&values.local_actual_bytes, &rank_actual, 1, MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Checkpoint actual bytes)");
  ++collectives;

  std::uint64_t manifest_size = 0U;
  std::uint64_t manifest_crc = 0U;
  std::uint64_t marker_size = 0U;
  bool manifest_ok = true;
  if (mpi.rank() == 0) {
    try {
      runtime::checkpoint_v2::Manifest manifest;
      manifest.rank_count = static_cast<std::uint32_t>(mpi.size());
      manifest.process_grid = decomposition.process_grid();
      manifest.fingerprints = identity;
      manifest.global_payload = global_payload;
      for (const auto &record : records) {
        manifest.ranks.push_back(
            {record.rank,
             {record.begin[0], record.begin[1], record.begin[2]},
             {record.end[0], record.end[1], record.end[2]},
             rank_filename(record.rank), record.logical, record.actual,
             record.crc, record.layout});
      }
      const auto bytes = runtime::checkpoint_v2::encode_manifest(manifest);
      const auto temp = directory / "manifest.v2.bin.tmp";
      const auto verified =
          runtime::checkpoint_v2::write_verified_temporary(temp, bytes);
      manifest_size = verified.actual_size;
      manifest_crc = verified.crc64;
      runtime::checkpoint_v2::publish_no_overwrite(
          temp, directory / "manifest.v2.bin");
    } catch (...) {
      manifest_ok = false;
    }
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest)");
  if (!status.ok) {
    values.manifest_crc = CheckpointV2CheckStatus::failed;
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::manifest, status.failing_rank);
  }
  runtime::check_mpi_result(
      MPI_Bcast(&manifest_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(Checkpoint manifest size)");
  ++collectives;
  runtime::check_mpi_result(
      MPI_Bcast(&manifest_crc, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(Checkpoint manifest CRC)");
  ++collectives;
  values.manifest_crc64 = manifest_crc;
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  ++values.crc_check_count;

  bool marker_ok = true;
  if (mpi.rank() == 0) {
    try {
      const runtime::checkpoint_v2::CompletedMarker marker{
          manifest_size, manifest_crc,
          common_fingerprint(identity, mpi.size(), decomposition.process_grid(),
                             global_payload)};
      const auto bytes =
          runtime::checkpoint_v2::encode_completed_marker(marker);
      marker_size = bytes.size();
      const auto temp = directory / "COMPLETED.tmp";
      static_cast<void>(
          runtime::checkpoint_v2::write_verified_temporary(temp, bytes));
      runtime::checkpoint_v2::publish_no_overwrite(temp,
                                                    directory / "COMPLETED");
    } catch (...) {
      marker_ok = false;
    }
  }
  status = converge(mpi, marker_ok, collectives,
                    "MPI_Allreduce(Checkpoint marker)");
  if (!status.ok) {
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::completed_marker, status.failing_rank);
  }
  runtime::check_mpi_result(
      MPI_Bcast(&marker_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(Checkpoint marker size)");
  ++collectives;
  ++values.crc_check_count;
  values.global_logical_bytes = global_logical;
  values.global_actual_bytes = rank_actual + manifest_size + marker_size;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.collective_count = collectives;
  values.publication = CheckpointV2CheckStatus::passed;
  values.disposition = CheckpointV2Disposition::completed;
  values.reason = CheckpointV2FailureReason::none;
  values.phase = CheckpointV2Phase::completed_marker;
  values.lowest_failing_rank = -1;
  return detail::CheckpointV2Access::make(values);
}

CheckpointV2ReadResult read_checkpoint_v2(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::FlowCaseConfig &config, FlowState &state,
    const std::filesystem::path &directory) {
  static_cast<void>(boundaries);
  detail::CheckpointV2ReportValues values;
  values.operation = CheckpointV2Operation::read;
  values.rank = mpi.rank();
  std::uint64_t collectives = 0U;
  bool entered = false;
  const auto fail = [&](CheckpointV2FailureReason reason,
                        CheckpointV2Phase phase, int rank) {
    values.disposition = CheckpointV2Disposition::failed;
    values.reason = reason;
    values.phase = phase;
    values.lowest_failing_rank = rank;
    values.collective_count = collectives;
    if (entered)
      values.rollback = CheckpointV2CheckStatus::passed;
    return detail::CheckpointV2Access::make_read(
        detail::CheckpointV2Access::make(values), {}, {}, false);
  };

  bool local_ok = true;
  bool common_path = false;
  try {
    common_path = path_agrees(mpi, directory, collectives);
  } catch (...) {
    throw;
  }
  try {
    local_ok = detail::FlowStateCheckpointAccess::live(state) &&
               !detail::FlowStateCheckpointAccess::attempt_active(state) &&
               detail::FlowStateCheckpointAccess::
                   diagnostic_identity_can_advance(state) &&
               geometry.compatible(topology) &&
               same(topology.global_extent(), decomposition.global_extent()) &&
               same(topology.owned_global_box(), decomposition.owned_box()) &&
               topology.local_face_count() ==
                   detail::FlowStateCheckpointAccess::layout(state).face_count &&
               same(decomposition.local_extent(),
                    detail::FlowStateCheckpointAccess::layout(state)
                        .cell_interior_extent) &&
               config.schema_version == 2 &&
               config.simulation_type ==
                   config::SimulationType::variable_density_flow;
    local_ok = local_ok && common_path;
  } catch (...) {
    local_ok = false;
  }
  auto status =
      converge(mpi, local_ok, collectives,
               "MPI_Allreduce(Checkpoint read preflight)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, status.failing_rank);

  const bool entry_ready =
      detail::FlowStateCheckpointAccess::read_transaction_ready(state);
  status = converge(mpi, entry_ready, collectives,
                    "MPI_Allreduce(Checkpoint read transaction readiness)");
  if (!status.ok) {
    values.transaction_entry = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::transaction_entry, status.failing_rank);
  }

  bool entered_local = true;
  try {
    detail::FlowStateCheckpointAccess::enter_read_transaction(state);
    entered = true;
  } catch (...) {
    entered_local = false;
  }
  status = converge(mpi, entered_local, collectives,
                    "MPI_Allreduce(Checkpoint read transaction entry)");
  if (!status.ok) {
    values.transaction_entry = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::transaction_entry, status.failing_rank);
  }
  values.transaction_entry = CheckpointV2CheckStatus::passed;

  runtime::checkpoint_v2::CompletedMarker marker;
  ByteVector manifest_bytes;
  runtime::checkpoint_v2::Manifest manifest;
  GlobalPayload global;
  bool marker_ok = true;
  try {
    const auto marker_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "COMPLETED", 40U);
    marker =
        runtime::checkpoint_v2::decode_completed_marker(marker_bytes);
  } catch (...) {
    marker_ok = false;
  }
  status = converge(mpi, marker_ok, collectives,
                    "MPI_Allreduce(Checkpoint marker read)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::file_integrity,
                CheckpointV2Phase::marker_read, status.failing_rank);

  bool manifest_ok = true;
  try {
    manifest_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "manifest.v2.bin", marker.manifest_actual_size);
    values.manifest_crc64 = runtime::checkpoint_v2::crc64_ecma(
        manifest_bytes.data(), manifest_bytes.size());
    if (values.manifest_crc64 != marker.manifest_crc64)
      throw runtime::Error("Checkpoint v2 manifest CRC is invalid");
    manifest = runtime::checkpoint_v2::decode_manifest(manifest_bytes);
    global = decode_global_payload(manifest.global_payload);
  } catch (...) {
    manifest_ok = false;
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest read)");
  if (!status.ok) {
    values.manifest_crc = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::file_integrity,
                CheckpointV2Phase::manifest_read, status.failing_rank);
  }
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  values.step = global.metadata.step;
  values.time_s = global.metadata.time_s;

  bool compatible = true;
  try {
    const auto &registry = detail::FlowStateCheckpointAccess::registry(state);
    const auto expected = fingerprints(decomposition, topology, config,
                                       registry, state.fields(), mpi.size());
    compatible =
        manifest.rank_count == static_cast<std::uint32_t>(mpi.size()) &&
        same(manifest.process_grid, decomposition.process_grid()) &&
        manifest.fingerprints == expected &&
        marker.common_fingerprint ==
            common_fingerprint(expected, mpi.size(),
                               decomposition.process_grid(),
                               manifest.global_payload) &&
        directory_inventory_valid(directory, mpi.size()) &&
        detail::TimeControlStateCodec::semantically_valid(
            config.time, config.density_model, global.metadata, global.time) &&
        (global.closure.has_value() ==
         (config.density_model == config::DensityModel::ideal_gas));
    const auto &record =
        manifest.ranks.at(static_cast<std::size_t>(mpi.rank()));
    compatible =
        compatible && record.rank == mpi.rank() &&
        same(runtime::Box3{record.owned_box_begin, record.owned_box_end},
             decomposition.owned_box()) &&
        record.local_layout_fingerprint ==
            local_layout_fingerprint(decomposition, topology, mpi.rank());
  } catch (...) {
    compatible = false;
  }
  status = converge(mpi, compatible, collectives,
                    "MPI_Allreduce(Checkpoint compatibility)");
  if (!status.ok) {
    values.fingerprint = CheckpointV2CheckStatus::failed;
    values.partition = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::layout,
                CheckpointV2Phase::manifest_read, status.failing_rank);
  }
  values.fingerprint = CheckpointV2CheckStatus::passed;
  values.partition = CheckpointV2CheckStatus::passed;

  ByteVector rank_bytes;
  std::pair<FlowLayerValues, FlowLayerValues> layers;
  bool rank_ok = true;
  try {
    const auto &record =
        manifest.ranks.at(static_cast<std::size_t>(mpi.rank()));
    rank_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / record.filename, record.actual_byte_size);
    values.local_actual_bytes = rank_bytes.size();
    values.local_crc64 = runtime::checkpoint_v2::crc64_ecma(
        rank_bytes.data(), rank_bytes.size());
    if (values.local_crc64 != record.crc64)
      throw runtime::Error("Checkpoint v2 rank CRC is invalid");
    const auto wrapper =
        runtime::checkpoint_v2::decode_rank_wrapper(rank_bytes);
    if (wrapper.rank != mpi.rank() || wrapper.rank_count != mpi.size())
      throw runtime::Error("Checkpoint v2 rank wrapper identity is invalid");
    layers = decode_rank_payload(wrapper.payload, state);
    values.local_logical_bytes = record.logical_byte_size;
  } catch (...) {
    rank_ok = false;
  }
  status = converge(mpi, rank_ok, collectives,
                    "MPI_Allreduce(Checkpoint rank read)");
  if (!status.ok) {
    values.rank_crc = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::file_integrity,
                CheckpointV2Phase::rank_read, status.failing_rank);
  }
  values.rank_crc = CheckpointV2CheckStatus::passed;

  bool state_ok =
      valid_layer(layers.first, config) &&
      valid_layer(layers.second, config);
  if (config.density_model == config::DensityModel::ideal_gas &&
      global.closure)
    state_ok =
        ideal_gas_layers_valid(mpi, topology, geometry, boundaries, config,
                               layers.first, layers.second, *global.closure,
                               collectives) &&
        state_ok;
  status = converge(mpi, state_ok, collectives,
                    "MPI_Allreduce(Checkpoint restored physical state)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare, status.failing_rank);

  std::optional<FlowState> replacement;
  bool prepared = true;
  try {
    replacement.emplace(
        detail::FlowStateCheckpointAccess::prepare_replacement(
            state, layers.first, layers.second, global.metadata));
  } catch (...) {
    prepared = false;
  }
  status = converge(mpi, prepared, collectives,
                    "MPI_Allreduce(Checkpoint restore prepare)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare, status.failing_rank);

  detail::FlowStateCheckpointAccess::publish_replacement(
      state, std::move(*replacement));
  values.publication = CheckpointV2CheckStatus::passed;
  values.rollback = CheckpointV2CheckStatus::not_checked;
  values.disposition = CheckpointV2Disposition::completed;
  values.reason = CheckpointV2FailureReason::none;
  values.phase = CheckpointV2Phase::restore_publish;
  values.lowest_failing_rank = -1;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
  std::uint64_t global_logical = 0U;
  std::uint64_t rank_actual = 0U;
  runtime::check_mpi_result(
      MPI_Allreduce(&values.local_logical_bytes, &global_logical, 1,
                    MPI_UINT64_T, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Checkpoint restored logical bytes)");
  ++collectives;
  runtime::check_mpi_result(
      MPI_Allreduce(&values.local_actual_bytes, &rank_actual, 1, MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Checkpoint restored actual bytes)");
  ++collectives;
  values.global_logical_bytes = global_logical;
  values.global_actual_bytes =
      rank_actual + marker.manifest_actual_size + 40U;
  values.collective_count = collectives;
  return detail::CheckpointV2Access::make_read(
      detail::CheckpointV2Access::make(values), global.time, global.closure,
      true);
}

CheckpointV2DiagnosticSource::CheckpointV2DiagnosticSource(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
CheckpointV2DiagnosticSource::~CheckpointV2DiagnosticSource() noexcept =
    default;
CheckpointV2DiagnosticSource::CheckpointV2DiagnosticSource(
    CheckpointV2DiagnosticSource &&) noexcept = default;
CheckpointV2DiagnosticSource
checkpoint_v2_diagnostic_source(const CheckpointV2Report &report) {
  return CheckpointV2DiagnosticSource(
      std::make_unique<CheckpointV2DiagnosticSource::Impl>(report));
}

} // namespace hundun::flow

// SPDX-License-Identifier: Apache-2.0

#include "checkpoint_v2_protocol.hpp"
#include "checkpoint_v2_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"
#include "hundun/flow/material_density_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/checkpoint_v2_product_oracle.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Bytes = std::vector<std::uint8_t>;

std::uint64_t bits_of(double value);

void append_u8(Bytes &bytes, std::uint8_t value) { bytes.push_back(value); }
void append_u32(Bytes &bytes, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U)
    append_u8(bytes, static_cast<std::uint8_t>(value >> shift));
}
void append_i32(Bytes &bytes, std::int32_t value) {
  std::uint32_t wire{};
  std::memcpy(&wire, &value, sizeof(wire));
  append_u32(bytes, wire);
}
void append_u64(Bytes &bytes, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    append_u8(bytes, static_cast<std::uint8_t>(value >> shift));
}
void append_f64(Bytes &bytes, double value) {
  std::uint64_t wire{};
  std::memcpy(&wire, &value, sizeof(wire));
  append_u64(bytes, wire);
}
void append_int3(Bytes &bytes, hundun::runtime::Int3 value) {
  append_i32(bytes, value.x);
  append_i32(bytes, value.y);
  append_i32(bytes, value.z);
}
void append_real3(Bytes &bytes, hundun::runtime::Real3 value) {
  append_f64(bytes, value.x);
  append_f64(bytes, value.y);
  append_f64(bytes, value.z);
}
void append_box3(Bytes &bytes, hundun::runtime::Box3 value) {
  append_int3(bytes, value.begin);
  append_int3(bytes, value.end);
}
void append_string(Bytes &bytes, const std::string &value) {
  append_u32(bytes, static_cast<std::uint32_t>(value.size()));
  bytes.insert(bytes.end(), value.begin(), value.end());
}
Bytes independent_global_payload(
    hundun::flow::AcceptedStepMetadata metadata,
    const hundun::flow::TimeControlState &time,
    const std::optional<hundun::flow::IdealGasClosureState> &closure) {
  Bytes bytes;
  append_u32(bytes, 1U);
  append_u64(bytes, metadata.step);
  append_f64(bytes, metadata.time_s);
  append_f64(bytes, metadata.dt_s);
  append_f64(bytes, metadata.previous_dt_s);
  append_u8(bytes,
            metadata.order == hundun::flow::MomentumTimeOrder::backward_euler
                ? 0U
                : 1U);
  append_u32(bytes, time.schema_version);
  append_u64(bytes, time.accepted_step);
  append_f64(bytes, time.proposed_next_dt_s);
  append_f64(bytes, time.last_accepted_dt_s);
  append_u8(bytes, time.last_accepted_order ==
                           hundun::flow::MomentumTimeOrder::backward_euler
                       ? 0U
                       : 1U);
  append_u8(bytes, time.history_ready ? 1U : 0U);
  append_u8(bytes, time.last_all_linear_solves_within_half_limit ? 1U : 0U);
  append_f64(bytes, time.last_convective_rate_per_s);
  append_f64(bytes, time.last_diffusive_rate_per_s);
  append_u8(bytes, time.last_stability_metrics_available ? 1U : 0U);
  append_u32(bytes, time.last_retry_count);
  append_u64(bytes, time.revision);
  append_u64(bytes, time.state_seal);
  append_u8(bytes, closure.has_value() ? 1U : 0U);
  if (closure) {
    append_u8(bytes, closure->mode ==
                             hundun::flow::IdealGasPressureMode::closed_dynamic
                         ? 0U
                         : 1U);
    append_f64(bytes, closure->thermodynamic_pressure_pa);
    append_u8(bytes, closure->target_mass_kg.has_value() ? 1U : 0U);
    if (closure->target_mass_kg)
      append_f64(bytes, *closure->target_mass_kg);
    append_u64(bytes, closure->revision);
  }
  return bytes;
}

void append_independent_rank_record(
    Bytes &bytes, std::uint8_t layer, std::uint8_t role,
    std::uint32_t transported, hundun::runtime::FieldId field,
    hundun::runtime::FunctionSpace space, std::uint32_t components,
    const std::vector<double> &values, std::uint64_t &logical) {
  append_u8(bytes, layer);
  append_u8(bytes, role);
  append_u32(bytes, transported);
  append_u32(bytes, field);
  append_u8(bytes,
            space == hundun::runtime::FunctionSpace::cell_average ? 0U : 1U);
  append_u8(bytes, 0U);
  append_u32(bytes, components);
  append_u64(bytes, values.size() / components);
  append_u64(bytes, values.size() * sizeof(double));
  for (const auto value : values)
    append_f64(bytes, value);
  logical += values.size() * sizeof(double);
}

Bytes independent_rank_payload(hundun::runtime::Int3 extent,
                               std::uint64_t face_count,
                               const hundun::flow::FlowFieldIds &fields,
                               const hundun::flow::FlowLayerValues &history,
                               const hundun::flow::FlowLayerValues &committed,
                               std::uint64_t &logical) {
  Bytes bytes;
  append_u32(bytes, 1U);
  append_int3(bytes, extent);
  append_u64(bytes, face_count);
  append_u32(bytes,
             static_cast<std::uint32_t>(fields.transported_cell_fields.size()));
  append_u32(bytes, static_cast<std::uint32_t>(
                        2U * (5U + fields.transported_cell_fields.size())));
  logical = 0U;
  const auto layer = [&](std::uint8_t layer_id,
                         const hundun::flow::FlowLayerValues &values) {
    append_independent_rank_record(bytes, layer_id, 0U, 0U, fields.density,
                                   hundun::runtime::FunctionSpace::cell_average,
                                   1U, values.density, logical);
    append_independent_rank_record(bytes, layer_id, 1U, 0U, fields.velocity,
                                   hundun::runtime::FunctionSpace::cell_average,
                                   3U, values.velocity, logical);
    append_independent_rank_record(bytes, layer_id, 2U, 0U,
                                   fields.mechanical_pressure,
                                   hundun::runtime::FunctionSpace::cell_average,
                                   1U, values.mechanical_pressure, logical);
    append_independent_rank_record(bytes, layer_id, 3U, 0U,
                                   fields.face_velocity,
                                   hundun::runtime::FunctionSpace::face_value,
                                   3U, values.face_velocity, logical);
    append_independent_rank_record(bytes, layer_id, 4U, 0U,
                                   fields.face_mass_flux,
                                   hundun::runtime::FunctionSpace::face_value,
                                   1U, values.face_mass_flux, logical);
    for (std::size_t index = 0; index < fields.transported_cell_fields.size();
         ++index)
      append_independent_rank_record(
          bytes, layer_id, 5U, static_cast<std::uint32_t>(index),
          fields.transported_cell_fields[index],
          hundun::runtime::FunctionSpace::cell_average, 1U,
          values.transported_cell_fields[index], logical);
  };
  layer(0U, history);
  layer(1U, committed);
  return bytes;
}
std::uint64_t independent_crc64(const Bytes &bytes) {
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  std::uint64_t crc{};
  for (const auto byte : bytes) {
    crc ^= static_cast<std::uint64_t>(byte) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
                                               : crc << 1U;
  }
  return crc;
}
template <class Authenticator>
void require_authenticated_mutation_matrix(const Bytes &bytes,
                                           Authenticator authenticate) {
  HUNDUN_CHECK(authenticate(bytes));
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    auto changed = bytes;
    changed[index] ^= 1U;
    HUNDUN_CHECK(!authenticate(changed));
  }
  for (std::size_t size = 0; size < bytes.size(); ++size) {
    auto changed = bytes;
    changed.resize(size);
    HUNDUN_CHECK(!authenticate(changed));
  }
  auto trailing = bytes;
  trailing.push_back(0U);
  HUNDUN_CHECK(!authenticate(trailing));
}
Bytes canonical_prefix(std::string_view domain) {
  Bytes bytes(domain.begin(), domain.end());
  append_u8(bytes, 0U);
  append_u32(bytes, 1U);
  return bytes;
}
void append_optional_double(Bytes &bytes, const std::optional<double> &value) {
  append_u8(bytes, value.has_value() ? 1U : 0U);
  if (value)
    append_f64(bytes, *value);
}
void append_optional_real3(Bytes &bytes,
                           const std::optional<hundun::runtime::Real3> &value) {
  append_u8(bytes, value.has_value() ? 1U : 0U);
  if (value)
    append_real3(bytes, *value);
}
std::uint64_t
independent_resolved_fingerprint(const hundun::config::FlowCaseConfig &config) {
  auto bytes = canonical_prefix("hundun.checkpoint-v2.resolved-case.v1");
  append_i32(bytes, config.schema_version);
  append_u8(bytes, static_cast<std::uint8_t>(config.simulation_type));
  append_u8(bytes, static_cast<std::uint8_t>(config.density_model));
  append_u8(bytes, config.resources.expected_ranks.has_value() ? 1U : 0U);
  if (config.resources.expected_ranks)
    append_i32(bytes, *config.resources.expected_ranks);
  append_u8(bytes, config.resources.process_grid.has_value() ? 1U : 0U);
  if (config.resources.process_grid)
    append_int3(bytes, *config.resources.process_grid);
  append_int3(bytes, config.mesh.cells);
  append_real3(bytes, config.mesh.origin_m);
  append_real3(bytes, config.mesh.length_m);
  append_u8(bytes, static_cast<std::uint8_t>(config.mesh.mapping));
  append_optional_real3(bytes, config.mesh.warp_amplitude);
  append_u8(bytes, static_cast<std::uint8_t>(config.time.mode));
  append_i32(bytes, config.time.steps);
  for (const auto value :
       {config.time.initial_dt_s, config.time.min_dt_s, config.time.max_dt_s,
        config.time.cfl_target, config.time.diffusion_number_target,
        config.time.growth_factor, config.time.retry_factor})
    append_f64(bytes, value);
  append_i32(bytes, config.time.max_retries);
  append_f64(bytes, config.physics.rho_ref_kg_per_m3);
  append_f64(bytes, config.physics.dynamic_viscosity_pa_s);
  append_f64(bytes, config.physics.inlet_consistency_rtol);
  append_optional_double(bytes, config.physics.cp_J_per_kg_K);
  append_optional_double(bytes, config.physics.gas_constant_J_per_kg_K);
  append_optional_double(bytes, config.physics.thermodynamic_pressure_pa);
  append_u64(bytes, config.scalars.size());
  for (const auto &scalar : config.scalars) {
    append_string(bytes, scalar.name);
    append_f64(bytes, scalar.diffusivity_m2_per_s);
  }
  for (const auto &patch : config.boundaries) {
    append_u8(bytes, static_cast<std::uint8_t>(patch.patch));
    append_u8(bytes, static_cast<std::uint8_t>(patch.type));
    append_optional_real3(bytes, patch.velocity_m_per_s);
    append_u8(bytes, patch.thermal_authority.has_value() ? 1U : 0U);
    if (patch.thermal_authority)
      append_u8(bytes, static_cast<std::uint8_t>(*patch.thermal_authority));
    append_optional_double(bytes, patch.temperature_K);
    append_optional_double(bytes, patch.enthalpy_J_per_kg);
    append_optional_double(bytes, patch.density_kg_per_m3);
    append_u8(bytes, patch.scalar_values.has_value() ? 1U : 0U);
    if (patch.scalar_values) {
      append_u64(bytes, patch.scalar_values->size());
      for (const auto &item : *patch.scalar_values) {
        append_string(bytes, item.name);
        append_f64(bytes, item.value);
      }
    }
    append_optional_double(bytes, patch.pressure_perturbation_pa);
  }
  return independent_crc64(bytes);
}

std::uint64_t independent_boundary_fingerprint(
    const hundun::boundary::BoundaryRegistry &boundaries) {
  auto bytes = canonical_prefix("hundun.checkpoint-v2.boundary.v1");
  append_u64(bytes, boundaries.scalar_count());
  for (std::size_t index = 0; index < boundaries.scalar_count(); ++index)
    append_string(bytes, std::string(boundaries.scalar_name(index)));
  append_u8(bytes, boundaries.open_domain() ? 1U : 0U);
  append_u8(bytes, boundaries.velocity_inlet_patch_id().has_value() ? 1U : 0U);
  if (boundaries.velocity_inlet_patch_id())
    append_u64(bytes, *boundaries.velocity_inlet_patch_id());
  append_u8(bytes, boundaries.pressure_outlet_patch_id().has_value() ? 1U : 0U);
  if (boundaries.pressure_outlet_patch_id())
    append_u64(bytes, *boundaries.pressure_outlet_patch_id());
  for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
    const auto &item = boundaries.patch(patch);
    append_u64(bytes, item.stable_id());
    append_string(bytes, std::string(item.name()));
    append_u8(bytes, static_cast<std::uint8_t>(item.kind()));
    append_u8(bytes, static_cast<std::uint8_t>(item.velocity_rule()));
    append_u8(bytes, static_cast<std::uint8_t>(item.pressure_rule()));
    append_u8(bytes, static_cast<std::uint8_t>(item.density_rule()));
    append_u8(bytes, static_cast<std::uint8_t>(item.enthalpy_rule()));
    append_u8(bytes, static_cast<std::uint8_t>(item.scalar_rule()));
    append_u8(bytes, static_cast<std::uint8_t>(item.mass_flux_rule()));
    append_u8(bytes, item.paired_patch_id().has_value() ? 1U : 0U);
    if (item.paired_patch_id())
      append_u64(bytes, *item.paired_patch_id());
    append_u8(bytes, item.inlet_state().has_value() ? 1U : 0U);
    if (item.inlet_state()) {
      append_real3(bytes, item.inlet_state()->velocity_m_per_s);
      append_f64(bytes, item.inlet_state()->density_kg_per_m3);
      append_f64(bytes, item.inlet_state()->enthalpy_J_per_kg);
      append_optional_double(bytes, item.inlet_state()->temperature_K);
      append_u64(bytes, item.inlet_state()->scalar_values.size());
      for (const auto value : item.inlet_state()->scalar_values)
        append_f64(bytes, value);
    }
    append_optional_double(bytes, item.pressure_value_pa());
  }
  return independent_crc64(bytes);
}

std::uint64_t independent_local_layout_fingerprint(
    const hundun::runtime::StructuredDecomposition &decomposition,
    const hundun::mesh::MeshTopology &topology, int rank, int rank_count) {
  auto bytes = canonical_prefix("hundun.checkpoint-v2.local-layout.v1");
  append_i32(bytes, rank);
  append_i32(bytes, rank_count);
  append_int3(bytes, decomposition.process_grid());
  append_box3(bytes, decomposition.owned_box());
  append_int3(bytes, decomposition.local_extent());
  append_u64(bytes, topology.local_cell_count());
  for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
    append_u8(bytes, topology.cell_ownership(cell) ==
                             hundun::mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
    append_u64(bytes, topology.global_cell_id(cell));
  }
  append_u64(bytes, topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    append_u8(bytes, topology.face_ownership(face) ==
                             hundun::mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
    append_u64(bytes, topology.global_face_id(face));
  }
  return independent_crc64(bytes);
}

std::uint64_t independent_local_topology_fingerprint(
    const hundun::runtime::StructuredDecomposition &decomposition,
    const hundun::mesh::MeshTopology &topology, int rank, int rank_count) {
  auto bytes = canonical_prefix("hundun.checkpoint-v2.topology-local.v1");
  append_i32(bytes, rank);
  append_i32(bytes, rank_count);
  append_int3(bytes, decomposition.process_grid());
  append_int3(bytes, topology.global_extent());
  append_box3(bytes, topology.owned_global_box());
  append_u64(bytes, topology.global_cell_count());
  append_u64(bytes, topology.global_face_count());
  for (const auto axis : {hundun::mesh::FaceAxis::x, hundun::mesh::FaceAxis::y,
                          hundun::mesh::FaceAxis::z})
    append_u64(bytes, topology.global_face_count(axis));
  append_u64(bytes, topology.owned_cell_count());
  append_u64(bytes, topology.ghost_cell_count());
  append_u64(bytes, topology.local_cell_count());
  append_u64(bytes, topology.owned_face_count());
  append_u64(bytes, topology.ghost_face_count());
  append_u64(bytes, topology.local_face_count());
  for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
    append_u64(bytes, topology.global_cell_id(cell));
    append_int3(bytes, topology.global_cell(cell));
    append_u8(bytes, topology.cell_ownership(cell) ==
                             hundun::mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
  }
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    append_u64(bytes, topology.global_face_id(face));
    const auto logical = topology.logical_face(face);
    append_u8(bytes, logical.axis == hundun::mesh::FaceAxis::x   ? 0U
                     : logical.axis == hundun::mesh::FaceAxis::y ? 1U
                                                                 : 2U);
    append_int3(bytes, logical.coordinate);
    append_u8(bytes, topology.face_ownership(face) ==
                             hundun::mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
    append_u64(bytes, topology.global_cell_id(topology.owner(face)));
    const auto neighbour = topology.neighbour(face);
    append_u8(bytes, neighbour.has_value() ? 1U : 0U);
    if (neighbour)
      append_u64(bytes, topology.global_cell_id(*neighbour));
    const auto patch = topology.patch_id(face);
    append_u8(bytes, patch.has_value() ? 1U : 0U);
    if (patch)
      append_u64(bytes, *patch);
    const auto pair = topology.periodic_pair(face);
    append_u8(bytes, pair.has_value() ? 1U : 0U);
    if (pair)
      append_u64(bytes, *pair);
  }
  for (const auto &patch : topology.patches()) {
    append_u64(bytes, patch.stable_id());
    append_string(bytes, std::string(patch.name()));
    append_u8(
        bytes,
        patch.pairing_kind() == hundun::mesh::PatchPairingKind::none ? 0U : 1U);
    append_u8(bytes, patch.paired_patch_id().has_value() ? 1U : 0U);
    if (patch.paired_patch_id())
      append_u64(bytes, *patch.paired_patch_id());
    append_u64(bytes, patch.local_faces().size());
    for (const auto face : patch.local_faces())
      append_u64(bytes, face);
  }
  return independent_crc64(bytes);
}

std::uint64_t independent_local_geometry_fingerprint(
    const hundun::mesh::MeshTopology &topology,
    const hundun::mesh::MeshGeometry &geometry,
    const hundun::config::FlowCaseConfig &config, int rank) {
  auto bytes = canonical_prefix("hundun.checkpoint-v2.geometry-local.v1");
  append_i32(bytes, rank);
  append_u8(bytes,
            geometry.mapping_kind() == hundun::mesh::MappingKind::uniform_box
                ? 0U
                : 1U);
  append_int3(bytes, geometry.global_extent());
  append_box3(bytes, geometry.owned_global_box());
  append_real3(bytes, geometry.origin_m());
  append_real3(bytes, geometry.length_m());
  const auto spacing = geometry.uniform_spacing_m();
  append_u8(bytes, spacing.has_value() ? 1U : 0U);
  if (spacing)
    append_real3(bytes, *spacing);
  append_u8(bytes, config.mesh.warp_amplitude.has_value() ? 1U : 0U);
  if (config.mesh.warp_amplitude)
    append_real3(bytes, *config.mesh.warp_amplitude);
  const auto box = topology.owned_global_box();
  const auto vertices =
      static_cast<std::uint64_t>(box.end.x - box.begin.x + 1) *
      static_cast<std::uint64_t>(box.end.y - box.begin.y + 1) *
      static_cast<std::uint64_t>(box.end.z - box.begin.z + 1);
  append_u64(bytes, vertices);
  for (int k = box.begin.z; k <= box.end.z; ++k)
    for (int j = box.begin.y; j <= box.end.y; ++j)
      for (int i = box.begin.x; i <= box.end.x; ++i) {
        const hundun::runtime::Int3 coordinate{i, j, k};
        append_int3(bytes, coordinate);
        append_real3(bytes, geometry.vertex_position_m(coordinate));
      }
  append_u64(bytes, topology.local_cell_count());
  for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
    append_u64(bytes, topology.global_cell_id(cell));
    append_real3(bytes, geometry.cell_center_m(cell));
    append_f64(bytes, geometry.cell_volume_m3(cell));
    append_f64(bytes, geometry.minimum_jacobian_determinant_m3(cell));
    if (topology.cell_ownership(cell) == hundun::mesh::EntityOwnership::owned)
      append_real3(bytes, geometry.cell_closure_m2(cell));
  }
  append_u64(bytes, topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    append_u64(bytes, topology.global_face_id(face));
    append_real3(bytes, geometry.face_center_m(face));
    append_real3(bytes, geometry.face_displacement_m(face));
    const auto owner =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    const auto neighbour =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::neighbour);
    append_real3(bytes, owner);
    append_real3(bytes, neighbour);
    HUNDUN_CHECK(bits_of(neighbour.x) == bits_of(-owner.x));
    HUNDUN_CHECK(bits_of(neighbour.y) == bits_of(-owner.y));
    HUNDUN_CHECK(bits_of(neighbour.z) == bits_of(-owner.z));
    append_f64(bytes, geometry.face_area_m2(face));
    append_f64(bytes, geometry.face_skewness(face));
    append_f64(bytes, geometry.face_non_orthogonality_degrees(face));
  }
  return independent_crc64(bytes);
}

std::pair<std::uint64_t, std::uint64_t> independent_common_mesh_fingerprints(
    const hundun::runtime::MpiContext &mpi,
    const hundun::runtime::StructuredDecomposition &decomposition,
    const hundun::mesh::MeshTopology &topology,
    const hundun::mesh::MeshGeometry &geometry,
    const hundun::config::FlowCaseConfig &config) {
  const auto local_topology = independent_local_topology_fingerprint(
      decomposition, topology, mpi.rank(), mpi.size());
  const auto local_geometry = independent_local_geometry_fingerprint(
      topology, geometry, config, mpi.rank());
  std::vector<std::uint64_t> topology_parts(
      static_cast<std::size_t>(mpi.size()));
  std::vector<std::uint64_t> geometry_parts(
      static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_topology, 1, MPI_UINT64_T,
                             topology_parts.data(), 1, MPI_UINT64_T,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allgather(&local_geometry, 1, MPI_UINT64_T,
                             geometry_parts.data(), 1, MPI_UINT64_T,
                             mpi.comm()) == MPI_SUCCESS);
  auto topology_bytes =
      canonical_prefix("hundun.checkpoint-v2.topology-common.v1");
  append_i32(topology_bytes, mpi.size());
  append_int3(topology_bytes, decomposition.process_grid());
  append_int3(topology_bytes, topology.global_extent());
  append_u64(topology_bytes, topology_parts.size());
  for (int rank = 0; rank < mpi.size(); ++rank) {
    append_i32(topology_bytes, rank);
    append_u64(topology_bytes, topology_parts[static_cast<std::size_t>(rank)]);
  }
  auto geometry_bytes =
      canonical_prefix("hundun.checkpoint-v2.geometry-common.v1");
  append_u8(geometry_bytes,
            geometry.mapping_kind() == hundun::mesh::MappingKind::uniform_box
                ? 0U
                : 1U);
  append_int3(geometry_bytes, geometry.global_extent());
  append_real3(geometry_bytes, geometry.origin_m());
  append_real3(geometry_bytes, geometry.length_m());
  const auto spacing = geometry.uniform_spacing_m();
  append_u8(geometry_bytes, spacing.has_value() ? 1U : 0U);
  if (spacing)
    append_real3(geometry_bytes, *spacing);
  append_u8(geometry_bytes, config.mesh.warp_amplitude.has_value() ? 1U : 0U);
  if (config.mesh.warp_amplitude)
    append_real3(geometry_bytes, *config.mesh.warp_amplitude);
  append_u64(geometry_bytes, geometry_parts.size());
  for (int rank = 0; rank < mpi.size(); ++rank) {
    append_i32(geometry_bytes, rank);
    append_u64(geometry_bytes, geometry_parts[static_cast<std::size_t>(rank)]);
  }
  return {independent_crc64(topology_bytes), independent_crc64(geometry_bytes)};
}
std::uint64_t read_u64_at(const Bytes &bytes, std::size_t offset) {
  std::uint64_t result{};
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    result |= static_cast<std::uint64_t>(bytes.at(offset + shift / 8U))
              << shift;
  return result;
}
std::uint64_t bits_of(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}
Bytes read_file_bytes(const std::filesystem::path &path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Task23 test could not read numerical file");
  return {std::istreambuf_iterator<char>(stream),
          std::istreambuf_iterator<char>()};
}
void write_file_bytes(const std::filesystem::path &path, const Bytes &bytes) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  if (!stream)
    throw std::runtime_error("Task23 test could not write numerical file");
  if (!bytes.empty())
    stream.write(reinterpret_cast<const char *>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
  stream.close();
  if (stream.fail())
    throw std::runtime_error("Task23 test numerical file write failed");
}

std::filesystem::path rank_file_path(const std::filesystem::path &directory,
                                     int rank) {
  const auto digits = std::to_string(rank);
  return directory /
         ("rank-" + std::string(6U - digits.size(), '0') + digits + ".v2.bin");
}

std::size_t rank_record_value_offset(std::size_t cells, std::size_t faces,
                                     std::size_t transported,
                                     std::size_t record_index,
                                     std::size_t value_index) {
  const auto records_per_layer = 5U + transported;
  if (record_index >= 2U * records_per_layer)
    throw std::runtime_error("Task23 test record index is invalid");
  const auto value_count = [&](std::size_t record) {
    switch (record % records_per_layer) {
    case 0U:
    case 2U:
      return cells;
    case 1U:
      return cells * 3U;
    case 3U:
      return faces * 3U;
    case 4U:
      return faces;
    default:
      return cells;
    }
  };
  std::size_t offset = 32U + 32U;
  for (std::size_t record = 0; record < record_index; ++record)
    offset += 32U + value_count(record) * sizeof(double);
  if (value_index >= value_count(record_index))
    throw std::runtime_error("Task23 test value index is invalid");
  return offset + 32U + value_index * sizeof(double);
}

void rebuild_checkpoint_rank_value(const hundun::runtime::MpiContext &mpi,
                                   const std::filesystem::path &directory,
                                   std::size_t cells, std::size_t faces,
                                   std::size_t transported,
                                   std::size_t record_index,
                                   std::uint64_t replacement_bits) {
  const int altered_rank = mpi.size() - 1;
  const auto own_rank_file = rank_file_path(directory, mpi.rank());
  if (mpi.rank() == altered_rank) {
    auto bytes = read_file_bytes(own_rank_file);
    const auto offset =
        rank_record_value_offset(cells, faces, transported, record_index, 0U);
    HUNDUN_CHECK(bytes.size() >= offset + sizeof(replacement_bits));
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes[offset + shift / 8U] =
          static_cast<std::uint8_t>(replacement_bits >> shift);
    write_file_bytes(own_rank_file, bytes);
  }
  mpi.barrier();
  const auto own_crc = independent_crc64(read_file_bytes(own_rank_file));
  std::vector<std::uint64_t> rank_crcs(static_cast<std::size_t>(mpi.size()));
  hundun::runtime::check_mpi_result(MPI_Allgather(&own_crc, 1, MPI_UINT64_T,
                                                  rank_crcs.data(), 1,
                                                  MPI_UINT64_T, mpi.comm()),
                                    "MPI_Allgather(Task23 rebuilt rank CRCs)");
  if (mpi.rank() == 0) {
    const auto manifest_path = directory / "manifest.v2.bin";
    auto manifest_bytes = read_file_bytes(manifest_path);
    const auto global_payload_size = manifest_bytes.size() - 84U -
                                     static_cast<std::size_t>(mpi.size()) * 82U;
    auto manifest = hundun::runtime::checkpoint_v2::decode_manifest(
        manifest_bytes, static_cast<std::uint32_t>(mpi.size()),
        global_payload_size);
    for (int rank = 0; rank < mpi.size(); ++rank)
      manifest.ranks[static_cast<std::size_t>(rank)].crc64 =
          rank_crcs[static_cast<std::size_t>(rank)];
    manifest_bytes = hundun::runtime::checkpoint_v2::encode_manifest(manifest);
    write_file_bytes(manifest_path, manifest_bytes);
    const auto marker_path = directory / "COMPLETED";
    auto marker = hundun::runtime::checkpoint_v2::decode_completed_marker(
        read_file_bytes(marker_path));
    marker.manifest_actual_size = manifest_bytes.size();
    marker.manifest_crc64 = independent_crc64(manifest_bytes);
    write_file_bytes(
        marker_path,
        hundun::runtime::checkpoint_v2::encode_completed_marker(marker));
  }
  mpi.barrier();
}

void rebuild_checkpoint_global_bytes(const hundun::runtime::MpiContext &mpi,
                                     const std::filesystem::path &directory,
                                     std::size_t global_payload_offset,
                                     const Bytes &replacement) {
  if (mpi.rank() == 0) {
    const auto manifest_path = directory / "manifest.v2.bin";
    auto manifest = read_file_bytes(manifest_path);
    HUNDUN_CHECK(manifest.size() >= 80U);
    const auto global_size = read_u64_at(manifest, 72U);
    HUNDUN_CHECK(global_payload_offset + replacement.size() <= global_size);
    const auto absolute = 80U + global_payload_offset;
    std::copy(replacement.begin(), replacement.end(),
              manifest.begin() + static_cast<std::ptrdiff_t>(absolute));
    write_file_bytes(manifest_path, manifest);

    Bytes common;
    for (std::size_t index = 0; index < 5U; ++index)
      append_u64(common, read_u64_at(manifest, 32U + index * 8U));
    std::uint32_t rank_bits{};
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      rank_bits |= static_cast<std::uint32_t>(manifest[16U + shift / 8U])
                   << shift;
    append_u32(common, rank_bits);
    common.insert(common.end(), manifest.begin() + 20, manifest.begin() + 32);
    common.insert(common.end(), manifest.begin() + 80,
                  manifest.begin() + 80 +
                      static_cast<std::ptrdiff_t>(global_size));

    Bytes marker{'H', 'F', 'C', '2', 'D', 'O', 'N', 0U};
    append_u32(marker, 2U);
    append_u32(marker, UINT32_C(0x01020304));
    append_u64(marker, manifest.size());
    append_u64(marker, independent_crc64(manifest));
    append_u64(marker, independent_crc64(common));
    write_file_bytes(directory / "COMPLETED", marker);
  }
  mpi.barrier();
}
void rebuild_checkpoint_global_u64(const hundun::runtime::MpiContext &mpi,
                                   const std::filesystem::path &directory,
                                   std::size_t global_payload_offset,
                                   std::uint64_t replacement) {
  Bytes bytes;
  append_u64(bytes, replacement);
  rebuild_checkpoint_global_bytes(mpi, directory, global_payload_offset, bytes);
}

std::uint64_t independent_field_schema_fingerprint(
    const hundun::runtime::FieldRegistry &registry,
    const hundun::flow::FlowFieldIds &fields) {
  Bytes bytes;
  constexpr std::string_view domain = "hundun.checkpoint-v2.field-schema.v1";
  bytes.insert(bytes.end(), domain.begin(), domain.end());
  append_u8(bytes, 0U);
  append_u32(bytes, 1U);
  std::vector<hundun::runtime::FieldId> ids{
      fields.density, fields.velocity, fields.mechanical_pressure,
      fields.face_velocity, fields.face_mass_flux};
  ids.insert(ids.end(), fields.transported_cell_fields.begin(),
             fields.transported_cell_fields.end());
  append_u64(bytes, ids.size());
  for (std::size_t index = 0; index < ids.size(); ++index) {
    const auto &descriptor = registry.descriptor(ids[index]);
    append_u8(bytes, index < 5U ? static_cast<std::uint8_t>(index) : 5U);
    append_u64(bytes, index < 5U ? 0U : index - 5U);
    append_u64(bytes, ids[index]);
    append_string(bytes, descriptor.name);
    append_string(bytes, descriptor.unit);
    append_string(bytes, descriptor.owner);
    append_u8(bytes, static_cast<std::uint8_t>(descriptor.space));
    append_u8(bytes, static_cast<std::uint8_t>(descriptor.scalar_type));
    append_u32(bytes, descriptor.components);
    append_i32(bytes, descriptor.ghost_width);
    append_u8(bytes, descriptor.conservative ? 1U : 0U);
    append_u8(bytes, static_cast<std::uint8_t>(descriptor.restart));
    append_u8(bytes, static_cast<std::uint8_t>(descriptor.output));
  }
  return independent_crc64(bytes);
}

struct ExpectedProductReport final {
  hundun::flow::detail::CheckpointV2ReportValues values;
  std::uint64_t semantic_fingerprint{};
};

std::uint64_t independent_report_fingerprint(
    const hundun::flow::detail::CheckpointV2ReportValues &v) {
  Bytes bytes;
  append_string(bytes, "hundun.checkpoint-v2.report-semantic.v1");
  append_u8(bytes, static_cast<std::uint8_t>(v.operation));
  append_u8(bytes, static_cast<std::uint8_t>(v.disposition));
  append_u8(bytes, static_cast<std::uint8_t>(v.reason));
  append_u8(bytes, static_cast<std::uint8_t>(v.phase));
  append_i32(bytes, v.rank);
  append_i32(bytes, v.lowest_failing_rank);
  append_u64(bytes, v.step);
  append_f64(bytes, v.time_s);
  append_u64(bytes, v.local_logical_bytes);
  append_u64(bytes, v.local_actual_bytes);
  append_u64(bytes, v.global_logical_bytes);
  append_u64(bytes, v.global_actual_bytes);
  append_u64(bytes, v.local_crc64);
  append_u64(bytes, v.manifest_crc64);
  append_u64(bytes, v.file_count);
  append_u64(bytes, v.crc_check_count);
  append_u64(bytes, v.collective_count);
  append_u8(bytes, static_cast<std::uint8_t>(v.rank_crc));
  append_u8(bytes, static_cast<std::uint8_t>(v.manifest_crc));
  append_u8(bytes, static_cast<std::uint8_t>(v.exact_size_eof));
  append_u8(bytes, static_cast<std::uint8_t>(v.fingerprint));
  append_u8(bytes, static_cast<std::uint8_t>(v.partition));
  append_u8(bytes, static_cast<std::uint8_t>(v.transaction_entry));
  append_u8(bytes, static_cast<std::uint8_t>(v.publication));
  append_u8(bytes, static_cast<std::uint8_t>(v.rollback));
  return independent_crc64(bytes);
}

hundun::flow::detail::CheckpointV2ReportValues observed_report_values(
    const hundun::flow::CheckpointV2Report &report) {
  hundun::flow::detail::CheckpointV2ReportValues v;
  v.operation = report.operation();
  v.disposition = report.disposition();
  v.reason = report.reason();
  v.phase = report.phase();
  v.rank = report.rank();
  v.lowest_failing_rank = report.lowest_failing_rank();
  v.step = report.step();
  v.time_s = report.time_s();
  v.local_logical_bytes = report.local_logical_bytes();
  v.local_actual_bytes = report.local_actual_bytes();
  v.global_logical_bytes = report.global_logical_bytes();
  v.global_actual_bytes = report.global_actual_bytes();
  v.local_crc64 = report.local_crc64();
  v.manifest_crc64 = report.manifest_crc64();
  v.file_count = report.file_count();
  v.crc_check_count = report.crc_check_count();
  v.collective_count = report.collective_count();
  v.rank_crc = report.rank_crc_status();
  v.manifest_crc = report.manifest_crc_status();
  v.exact_size_eof = report.exact_size_and_eof_status();
  v.fingerprint = report.fingerprint_status();
  v.partition = report.partition_status();
  v.transaction_entry = report.transaction_entry_status();
  v.publication = report.publication_status();
  v.rollback = report.rollback_status();
  return v;
}

bool product_report_matches(
    const hundun::flow::CheckpointV2Report &report,
    const ExpectedProductReport &expected) {
  const auto &v = expected.values;
  return report.operation() == v.operation &&
         report.disposition() == v.disposition &&
         report.reason() == v.reason && report.phase() == v.phase &&
         report.rank() == v.rank &&
         report.lowest_failing_rank() == v.lowest_failing_rank &&
         report.step() == v.step && bits_of(report.time_s()) == bits_of(v.time_s) &&
         report.local_logical_bytes() == v.local_logical_bytes &&
         report.local_actual_bytes() == v.local_actual_bytes &&
         report.global_logical_bytes() == v.global_logical_bytes &&
         report.global_actual_bytes() == v.global_actual_bytes &&
         report.local_crc64() == v.local_crc64 &&
         report.manifest_crc64() == v.manifest_crc64 &&
         report.file_count() == v.file_count &&
         report.crc_check_count() == v.crc_check_count &&
         report.collective_count() == v.collective_count &&
         report.rank_crc_status() == v.rank_crc &&
         report.manifest_crc_status() == v.manifest_crc &&
         report.exact_size_and_eof_status() == v.exact_size_eof &&
         report.fingerprint_status() == v.fingerprint &&
         report.partition_status() == v.partition &&
         report.transaction_entry_status() == v.transaction_entry &&
         report.publication_status() == v.publication &&
         report.rollback_status() == v.rollback &&
         report.semantic_fingerprint() == expected.semantic_fingerprint;
}

void require_exact_product_report(
    const hundun::flow::CheckpointV2Report &report,
    const ExpectedProductReport &expected) {
  if (report.semantic_fingerprint() != expected.semantic_fingerprint)
    throw std::runtime_error("Task23 exact report semantic mismatch: actual=" +
                             std::to_string(report.semantic_fingerprint()) +
                             " expected=" +
                             std::to_string(expected.semantic_fingerprint) +
                             " collective=" +
                             std::to_string(report.collective_count()) + "/" +
                             std::to_string(
                                 expected.values.collective_count) +
                             " observed-recomputed=" +
                             std::to_string(independent_report_fingerprint(
                                 observed_report_values(report))));
  if (report.global_actual_bytes() != expected.values.global_actual_bytes)
    throw std::runtime_error("Task23 exact report global actual mismatch");
  if (report.collective_count() != expected.values.collective_count)
    throw std::runtime_error("Task23 exact report collective count mismatch");
  HUNDUN_CHECK(product_report_matches(report, expected));
  auto changed_count = expected;
  ++changed_count.values.collective_count;
  HUNDUN_CHECK(!product_report_matches(report, changed_count));
  auto changed_tri_state = expected;
  changed_tri_state.values.manifest_crc =
      changed_tri_state.values.manifest_crc ==
              hundun::flow::CheckpointV2CheckStatus::failed
          ? hundun::flow::CheckpointV2CheckStatus::passed
          : hundun::flow::CheckpointV2CheckStatus::failed;
  HUNDUN_CHECK(!product_report_matches(report, changed_tri_state));
  auto changed_identity = expected;
  changed_identity.semantic_fingerprint ^= UINT64_C(1);
  HUNDUN_CHECK(!product_report_matches(report, changed_identity));
}

void require_consistent_product_report(const hundun::runtime::MpiContext &mpi,
                           const hundun::flow::CheckpointV2Report &report) {
  std::uint64_t time_bits{};
  const double time_s = report.time_s();
  std::memcpy(&time_bits, &time_s, sizeof(time_bits));
  const std::array<std::uint64_t, 19> local{
      static_cast<std::uint8_t>(report.operation()),
      static_cast<std::uint8_t>(report.disposition()),
      static_cast<std::uint8_t>(report.reason()),
      static_cast<std::uint8_t>(report.phase()),
      static_cast<std::uint64_t>(report.lowest_failing_rank()),
      report.step(),
      time_bits,
      report.global_logical_bytes(),
      report.global_actual_bytes(),
      report.manifest_crc64(),
      report.file_count(),
      report.crc_check_count(),
      report.collective_count(),
      static_cast<std::uint8_t>(report.manifest_crc_status()),
      static_cast<std::uint8_t>(report.exact_size_and_eof_status()),
      static_cast<std::uint8_t>(report.fingerprint_status()),
      static_cast<std::uint8_t>(report.partition_status()),
      static_cast<std::uint8_t>(report.transaction_entry_status()),
      (static_cast<std::uint64_t>(report.publication_status()) << 8U) |
          static_cast<std::uint8_t>(report.rollback_status())};
  auto authority = local;
  hundun::runtime::check_mpi_result(
      MPI_Bcast(authority.data(), static_cast<int>(authority.size()),
                MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(Task23 common report oracle)");
  HUNDUN_CHECK(local == authority);

  class ProductReportSink final : public hundun::diagnostics::DiagnosticSink {
  public:
    void
    submit(const hundun::diagnostics::DiagnosticRecord &candidate) override {
      record = candidate;
      ++calls;
    }
    hundun::diagnostics::DiagnosticRecord record;
    int calls{};
  };
  const auto phase = [&] {
    using Phase = hundun::flow::CheckpointV2Phase;
    switch (report.phase()) {
    case Phase::none:
      return "none";
    case Phase::preflight:
      return "preflight";
    case Phase::transaction_entry:
      return "transaction-entry";
    case Phase::rank_payload:
      return "rank-payload";
    case Phase::rank_temporary_file:
      return "rank-temporary-file";
    case Phase::rank_publish:
      return "rank-publish";
    case Phase::manifest:
      return "manifest";
    case Phase::completed_marker:
      return "completed-marker";
    case Phase::marker_read:
      return "marker-read";
    case Phase::manifest_read:
      return "manifest-read";
    case Phase::rank_read:
      return "rank-read";
    case Phase::restore_prepare:
      return "restore-prepare";
    case Phase::restore_publish:
      return "restore-publish";
    }
    return "none";
  }();
  const auto source = hundun::flow::checkpoint_v2_diagnostic_source(report);
  for (const auto level :
       {hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticLevel::invariants,
        hundun::diagnostics::DiagnosticLevel::counters}) {
    const hundun::diagnostics::DiagnosticRequest diagnostic_request{
        level,
        hundun::diagnostics::DiagnosticScope::local,
        {report.rank(), report.step(), report.time_s(), phase},
        {},
        0U};
    ProductReportSink diagnostic;
    hundun::diagnostics::collect_diagnostics(source, diagnostic_request,
                                             diagnostic);
    HUNDUN_CHECK(diagnostic.calls == 1);
    hundun::test::checkpoint_v2_oracle::
        require_exact_product_diagnostic_record(report, level,
                                                diagnostic.record);
  }
}

void require_destination_report_authority(
    const hundun::runtime::MpiContext &mpi,
    const hundun::flow::test::CheckpointV2DeepSnapshot &before,
    const hundun::flow::CheckpointV2Report &report) {
  HUNDUN_CHECK(report.step() == before.metadata.step);
  HUNDUN_CHECK(bits_of(report.time_s()) == bits_of(before.metadata.time_s));
  require_consistent_product_report(mpi, report);
  HUNDUN_CHECK(report.semantic_fingerprint() ==
               independent_report_fingerprint(observed_report_values(report)));
}

ExpectedProductReport expected_read_preflight_failure_report(
    const hundun::runtime::MpiContext &mpi,
    const hundun::flow::test::CheckpointV2DeepSnapshot &before) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation = hundun::flow::CheckpointV2Operation::read;
  values.disposition = hundun::flow::CheckpointV2Disposition::failed;
  values.reason = hundun::flow::CheckpointV2FailureReason::invalid_input;
  values.phase = hundun::flow::CheckpointV2Phase::preflight;
  values.rank = mpi.rank();
  values.lowest_failing_rank = 0;
  values.step = before.metadata.step;
  values.time_s = before.metadata.time_s;
  values.local_logical_bytes = 0U;
  values.local_actual_bytes = 0U;
  values.global_logical_bytes = 0U;
  values.global_actual_bytes = 0U;
  values.local_crc64 = 0U;
  values.manifest_crc64 = 0U;
  values.file_count = 0U;
  values.crc_check_count = 0U;
  values.collective_count = 7U;
  values.rank_crc = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.partition = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.transaction_entry = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.publication = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.rollback = hundun::flow::CheckpointV2CheckStatus::not_checked;
  return {values, independent_report_fingerprint(values)};
}

ExpectedProductReport expected_constructible_fingerprint_failure_report(
    const hundun::runtime::MpiContext &mpi,
    const std::filesystem::path &checkpoint_directory,
    const hundun::flow::test::CheckpointV2DeepSnapshot &before) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation = hundun::flow::CheckpointV2Operation::read;
  values.disposition = hundun::flow::CheckpointV2Disposition::failed;
  values.reason =
      hundun::flow::CheckpointV2FailureReason::file_integrity;
  values.phase = hundun::flow::CheckpointV2Phase::manifest_read;
  values.rank = mpi.rank();
  values.lowest_failing_rank = 0;
  values.step = before.metadata.step;
  values.time_s = before.metadata.time_s;
  values.local_logical_bytes = 0U;
  values.local_actual_bytes = 0U;
  values.global_logical_bytes = 0U;
  values.global_actual_bytes = 0U;
  values.local_crc64 = 0U;
  values.manifest_crc64 = independent_crc64(
      read_file_bytes(checkpoint_directory / "manifest.v2.bin"));
  values.file_count = 2U;
  values.crc_check_count = 1U;
  values.collective_count = 20U;
  values.rank_crc =
      hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::passed;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::failed;
  values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
  values.transaction_entry =
      hundun::flow::CheckpointV2CheckStatus::passed;
  values.publication =
      hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.rollback = hundun::flow::CheckpointV2CheckStatus::passed;
  return {values, independent_report_fingerprint(values)};
}

void require_constructible_fingerprint_failure(
    const hundun::runtime::MpiContext &mpi,
    const std::filesystem::path &checkpoint_directory,
    const hundun::flow::test::CheckpointV2DeepSnapshot &before,
    const hundun::flow::FlowState &state,
    const hundun::flow::CheckpointV2ReadResult &result) {
  HUNDUN_CHECK(!result.restored());
  const auto expected = expected_constructible_fingerprint_failure_report(
      mpi, checkpoint_directory, before);
  require_consistent_product_report(mpi, result.report());
  require_exact_product_report(result.report(), expected);
  HUNDUN_CHECK(
      hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
          before,
          hundun::flow::test::CheckpointV2TestAccess::snapshot(state)));
}

template <class Function> bool rejects(Function &&function) {
  try {
    function();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

hundun::config::FlowCaseConfig make_case(int ranks) {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "checkpoint-v2";
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::material;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid(ranks);
  config.mesh.cells = {8, 4, 3};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = hundun::config::MeshMapping::uniform_box;
  config.time = {hundun::config::TimeMode::fixed,
                 4,
                 0.01,
                 0.001,
                 0.1,
                 0.5,
                 0.25,
                 1.25,
                 0.5,
                 8};
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars = {{"alpha", 0.0}};
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, std::uint32_t n,
                                      bool conservative) {
  return {name,
          "1",
          "checkpoint-v2",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          n,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor face(const char *name, std::uint32_t n) {
  return {name,
          "1",
          "checkpoint-v2",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          n,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor
physical_cell(const char *name, const char *unit, bool conservative) {
  auto result = cell(name, 1U, conservative);
  result.unit = unit;
  return result;
}

void run(const hundun::runtime::MpiContext &mpi, bool acceptance) {
  using CheckpointAccess = hundun::flow::test::CheckpointV2TestAccess;
  HUNDUN_CHECK(hundun::test::
                   checkpoint_v2_state_equality_oracle_is_mutation_sensitive());
  auto config = make_case(mpi.size());
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, config.mesh.cells, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology, hundun::mesh::UniformBoxMapping(config.mesh.origin_m,
                                                config.mesh.length_m));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(config, topology);
  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(physical_cell("rho", "kg/m3", true));
  fields.velocity = registry.declare_field(cell("u", 3U, false));
  fields.mechanical_pressure = registry.declare_field(cell("pi", 1U, false));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields = {
      registry.declare_field(physical_cell("rho_h", "J/m3", true)),
      registry.declare_field(physical_cell("rho_alpha", "kg/m3", true))};
  registry.freeze();
  const auto metadata = hundun::flow::AcceptedStepMetadata{
      0U, 0.0, config.time.initial_dt_s, 0.0,
      hundun::flow::MomentumTimeOrder::backward_euler};
  const auto destination_metadata = hundun::flow::AcceptedStepMetadata{
      17U, 1.25, config.time.initial_dt_s, 0.0,
      hundun::flow::MomentumTimeOrder::backward_euler};
  auto make_state = [&] {
    return hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields, metadata);
  };
  auto make_destination_state = [&] {
    return hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields, destination_metadata);
  };
  const std::size_t cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues history;
  history.density.resize(cells);
  history.velocity.resize(cells * 3U);
  history.mechanical_pressure.resize(cells);
  history.face_velocity.resize(topology.local_face_count() * 3U);
  history.face_mass_flux.resize(topology.local_face_count());
  history.transported_cell_fields = {std::vector<double>(cells),
                                     std::vector<double>(cells)};
  for (std::size_t cell_id = 0; cell_id < cells; ++cell_id) {
    history.density[cell_id] = 1.0 + 0.001 * (static_cast<double>(mpi.rank()) +
                   static_cast<double>(cell_id));
    history.mechanical_pressure[cell_id] = -0.25 + static_cast<double>(cell_id);
    history.transported_cell_fields[0][cell_id] =
        history.density[cell_id] * 300000.0;
    history.transported_cell_fields[1][cell_id] =
        history.density[cell_id] * 0.2;
    for (std::size_t component = 0; component < 3U; ++component)
      history.velocity[cell_id * 3U + component] =
          10.0 * static_cast<double>(mpi.rank()) +
          3.0 * static_cast<double>(cell_id) + static_cast<double>(component);
  }
  for (std::size_t face_id = 0; face_id < topology.local_face_count();
       ++face_id) {
    history.face_mass_flux[face_id] =
        (static_cast<double>(mpi.rank()) + static_cast<double>(face_id)) * 0.01;
    for (std::size_t component = 0; component < 3U; ++component)
      history.face_velocity[face_id * 3U + component] =
          static_cast<double>(face_id + component) * 0.02;
  }
  history.mechanical_pressure.front() = -0.0;
  auto committed = history;
  for (double &item : committed.density)
    item += 0.25;
  for (double &item : committed.transported_cell_fields[0])
    item += 0.125;

  auto source = make_state();
  source.seed_accepted_layers(history, committed);
  auto controller = hundun::flow::Bdf2RetryController::create(
      config.time, config.density_model, topology, geometry, mpi, source);
  const auto controller_state = controller.state();

  const auto directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-checkpoint-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(directory);
  mpi.barrier();
  const auto written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller_state, std::nullopt, directory);
  if (written.disposition() != hundun::flow::CheckpointV2Disposition::completed)
    throw std::runtime_error(
        "Checkpoint write failed: reason=" +
        std::to_string(static_cast<unsigned>(written.reason())) +
        " phase=" + std::to_string(static_cast<unsigned>(written.phase())) +
        " rank=" + std::to_string(written.lowest_failing_rank()));
  HUNDUN_CHECK(written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  HUNDUN_CHECK(written.file_count() ==
               static_cast<std::uint64_t>(mpi.size()) + 2U);
  HUNDUN_CHECK(std::filesystem::is_regular_file(directory / "COMPLETED"));
  const auto independent_global =
      independent_global_payload(metadata, controller_state, std::nullopt);
  const auto product_global =
      hundun::flow::test::checkpoint_v2_encode_global_payload_for_test(
          metadata, controller_state, std::nullopt);
  HUNDUN_CHECK(product_global == independent_global);
  std::uint64_t independent_rank_logical{};
  const auto independent_rank = independent_rank_payload(
      decomposition.local_extent(), topology.local_face_count(), fields,
      history, committed, independent_rank_logical);
  std::uint64_t product_rank_logical{};
  const auto product_rank =
      hundun::flow::test::checkpoint_v2_encode_rank_payload_for_test(
          source, product_rank_logical);
  HUNDUN_CHECK(product_rank == independent_rank);
  HUNDUN_CHECK(product_rank_logical == independent_rank_logical);
  const auto actual_rank_wrapper =
      read_file_bytes(rank_file_path(directory, mpi.rank()));
  HUNDUN_CHECK(actual_rank_wrapper.size() == independent_rank.size() + 32U);
  HUNDUN_CHECK(std::equal(independent_rank.begin(), independent_rank.end(),
                          actual_rank_wrapper.begin() + 32));
  std::uint64_t expected_global_logical{};
  std::uint64_t expected_rank_actual = independent_rank.size() + 32U;
  std::uint64_t expected_rank_actual_sum{};
  HUNDUN_CHECK(MPI_Allreduce(&independent_rank_logical,
                             &expected_global_logical, 1, MPI_UINT64_T,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&expected_rank_actual,
                             &expected_rank_actual_sum, 1, MPI_UINT64_T,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  const auto manifest_bytes_for_report =
      read_file_bytes(directory / "manifest.v2.bin");
  const auto expected_manifest_size =
      84U + independent_global.size() +
      82U * static_cast<std::size_t>(mpi.size());
  HUNDUN_CHECK(manifest_bytes_for_report.size() == expected_manifest_size);
  hundun::flow::detail::CheckpointV2ReportValues expected_write;
  expected_write.operation = hundun::flow::CheckpointV2Operation::write;
  expected_write.disposition =
      hundun::flow::CheckpointV2Disposition::completed;
  expected_write.reason = hundun::flow::CheckpointV2FailureReason::none;
  expected_write.phase = hundun::flow::CheckpointV2Phase::completed_marker;
  expected_write.rank = mpi.rank();
  expected_write.lowest_failing_rank = -1;
  expected_write.step = metadata.step;
  expected_write.time_s = metadata.time_s;
  expected_write.local_logical_bytes = independent_rank_logical;
  expected_write.local_actual_bytes = expected_rank_actual;
  expected_write.global_logical_bytes = expected_global_logical;
  expected_write.global_actual_bytes =
      expected_rank_actual_sum + expected_manifest_size + 40U;
  expected_write.local_crc64 = independent_crc64(actual_rank_wrapper);
  expected_write.manifest_crc64 =
      independent_crc64(manifest_bytes_for_report);
  expected_write.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  expected_write.crc_check_count =
      static_cast<std::uint64_t>(mpi.size()) + 2U;
  expected_write.collective_count = 32U;
  expected_write.rank_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.manifest_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.exact_size_eof =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.fingerprint =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.partition =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.transaction_entry =
      hundun::flow::CheckpointV2CheckStatus::not_checked;
  expected_write.publication =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_write.rollback =
      hundun::flow::CheckpointV2CheckStatus::not_checked;
  require_exact_product_report(
      written,
      {expected_write, independent_report_fingerprint(expected_write)});
  const auto independent_mesh_fingerprints =
      independent_common_mesh_fingerprints(mpi, decomposition, topology,
                                           geometry, config);
  const auto independent_layout = independent_local_layout_fingerprint(
      decomposition, topology, mpi.rank(), mpi.size());
  std::vector<std::uint64_t> independent_layouts(
      static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&independent_layout, 1, MPI_UINT64_T,
                             independent_layouts.data(), 1, MPI_UINT64_T,
                             mpi.comm()) == MPI_SUCCESS);
  if (mpi.rank() == 0) {
    std::ifstream manifest_stream(directory / "manifest.v2.bin",
                                  std::ios::binary);
    const Bytes manifest_bytes{std::istreambuf_iterator<char>(manifest_stream),
                               std::istreambuf_iterator<char>()};
    HUNDUN_CHECK(manifest_bytes.size() >= 72U);
    HUNDUN_CHECK(manifest_bytes.size() >= 80U + independent_global.size());
    HUNDUN_CHECK(std::equal(independent_global.begin(),
                            independent_global.end(),
                            manifest_bytes.begin() + 80));
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 32U) ==
                 independent_resolved_fingerprint(config));
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 40U) ==
                 independent_mesh_fingerprints.first);
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 48U) ==
                 independent_mesh_fingerprints.second);
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 56U) ==
                 independent_boundary_fingerprint(boundaries));
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 64U) ==
                 independent_field_schema_fingerprint(registry, fields));
    const auto first_record = 84U + independent_global.size();
    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto record = first_record + static_cast<std::size_t>(rank) * 82U;
      HUNDUN_CHECK(read_u64_at(manifest_bytes, record + 74U) ==
                   independent_layouts[static_cast<std::size_t>(rank)]);
    }
  }
  const auto baseline_schema_fingerprint =
      independent_field_schema_fingerprint(registry, fields);
  using DescriptorMutation =
      std::function<void(hundun::runtime::FieldDescriptor &)>;
  const std::vector<DescriptorMutation> descriptor_mutations{
      [](auto &item) { item.name += "-changed"; },
      [](auto &item) { item.unit += "-changed"; },
      [](auto &item) { item.owner += "-changed"; },
      [](auto &item) {
        item.space = hundun::runtime::FunctionSpace::vertex_value;
      },
      [](auto &item) { item.scalar_type = hundun::runtime::ScalarType::int32; },
      [](auto &item) { item.components += 1U; },
      [](auto &item) { item.ghost_width += 1; },
      [](auto &item) { item.conservative = !item.conservative; },
      [](auto &item) {
        item.restart = hundun::runtime::RestartPolicy::transient;
      },
      [](auto &item) {
        item.output = hundun::runtime::OutputPolicy::selected;
      }};
  const std::array ordered_schema_ids{fields.density,
                                      fields.velocity,
                                      fields.mechanical_pressure,
                                      fields.face_velocity,
                                      fields.face_mass_flux,
                                      fields.transported_cell_fields[0],
                                      fields.transported_cell_fields[1]};
  for (std::size_t descriptor_mutation_index = 0U;
       descriptor_mutation_index < descriptor_mutations.size();
       ++descriptor_mutation_index) {
    const auto &mutate = descriptor_mutations[descriptor_mutation_index];
    for (std::size_t changed_index = 0;
         changed_index < ordered_schema_ids.size(); ++changed_index) {
      std::array<hundun::runtime::FieldDescriptor, 7> descriptors;
      for (std::size_t index = 0; index < descriptors.size(); ++index)
        descriptors[index] = registry.descriptor(ordered_schema_ids[index]);
      mutate(descriptors[changed_index]);
      hundun::runtime::FieldRegistry changed_registry;
      hundun::flow::FlowFieldIds changed_fields;
      changed_fields.density = changed_registry.declare_field(descriptors[0]);
      changed_fields.velocity = changed_registry.declare_field(descriptors[1]);
      changed_fields.mechanical_pressure =
          changed_registry.declare_field(descriptors[2]);
      changed_fields.face_velocity =
          changed_registry.declare_field(descriptors[3]);
      changed_fields.face_mass_flux =
          changed_registry.declare_field(descriptors[4]);
      changed_fields.transported_cell_fields = {
          changed_registry.declare_field(descriptors[5]),
          changed_registry.declare_field(descriptors[6])};
      changed_registry.freeze();
      const auto independent_changed_schema =
          independent_field_schema_fingerprint(changed_registry,
                                               changed_fields);
      HUNDUN_CHECK(independent_changed_schema != baseline_schema_fingerprint);
      HUNDUN_CHECK(
          hundun::flow::test::checkpoint_v2_field_schema_fingerprint_for_test(
              changed_registry, changed_fields) == independent_changed_schema);
      if (changed_index != 4U &&
          (descriptor_mutation_index <= 2U ||
           descriptor_mutation_index >= 7U)) {
        auto changed_state = hundun::flow::FlowState::create(
            changed_registry,
            {decomposition.local_extent(), topology.local_face_count()},
            changed_fields, destination_metadata);
        changed_state.seed_accepted_layers(history, committed);
        const auto changed_before = CheckpointAccess::snapshot(changed_state);
        const auto changed_read = hundun::flow::read_checkpoint_v2(
            mpi, decomposition, topology, geometry, boundaries, config,
            changed_state, directory);
        require_constructible_fingerprint_failure(
            mpi, directory, changed_before, changed_state, changed_read);
      }
    }
  }
  {
    hundun::runtime::FieldRegistry shifted_registry;
    static_cast<void>(shifted_registry.declare_field(
        cell("checkpoint-fingerprint-id-shift", 1U, false)));
    hundun::flow::FlowFieldIds shifted_fields;
    shifted_fields.density =
        shifted_registry.declare_field(registry.descriptor(fields.density));
    shifted_fields.velocity =
        shifted_registry.declare_field(registry.descriptor(fields.velocity));
    shifted_fields.mechanical_pressure = shifted_registry.declare_field(
        registry.descriptor(fields.mechanical_pressure));
    shifted_fields.face_velocity = shifted_registry.declare_field(
        registry.descriptor(fields.face_velocity));
    shifted_fields.face_mass_flux = shifted_registry.declare_field(
        registry.descriptor(fields.face_mass_flux));
    for (const auto field : fields.transported_cell_fields)
      shifted_fields.transported_cell_fields.push_back(
          shifted_registry.declare_field(registry.descriptor(field)));
    shifted_registry.freeze();
    HUNDUN_CHECK(
        independent_field_schema_fingerprint(shifted_registry,
                                             shifted_fields) !=
        baseline_schema_fingerprint);
    auto shifted_state = hundun::flow::FlowState::create(
        shifted_registry,
        {decomposition.local_extent(), topology.local_face_count()},
        shifted_fields, destination_metadata);
    shifted_state.seed_accepted_layers(history, committed);
    const auto shifted_before = CheckpointAccess::snapshot(shifted_state);
    const auto shifted_read = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config,
        shifted_state, directory);
    require_constructible_fingerprint_failure(
        mpi, directory, shifted_before, shifted_state, shifted_read);
  }

  auto destination = make_destination_state();
  auto different = history;
  for (double &item : different.density)
    item += 9.0;
  destination.seed_accepted_layers(different, different);
  const auto destination_before = CheckpointAccess::snapshot(destination);
  const auto restored = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, destination,
      directory);
  HUNDUN_CHECK(restored.restored());
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      history, destination.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      committed, destination.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      committed, destination.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      metadata, destination.metadata()));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      controller_state, restored.time_control_state()));
  HUNDUN_CHECK(
      CheckpointAccess::all_cell_ghosts_are_positive_zero(destination));
  require_destination_report_authority(mpi, destination_before,
                                       restored.report());
  hundun::flow::detail::CheckpointV2ReportValues expected_read;
  expected_read.operation = hundun::flow::CheckpointV2Operation::read;
  expected_read.disposition =
      hundun::flow::CheckpointV2Disposition::completed;
  expected_read.reason = hundun::flow::CheckpointV2FailureReason::none;
  expected_read.phase = hundun::flow::CheckpointV2Phase::restore_publish;
  expected_read.rank = mpi.rank();
  expected_read.lowest_failing_rank = -1;
  expected_read.step = destination_metadata.step;
  expected_read.time_s = destination_metadata.time_s;
  expected_read.local_logical_bytes = independent_rank_logical;
  expected_read.local_actual_bytes = expected_rank_actual;
  expected_read.global_logical_bytes = expected_global_logical;
  expected_read.global_actual_bytes =
      expected_rank_actual_sum + expected_manifest_size + 40U;
  expected_read.local_crc64 = independent_crc64(actual_rank_wrapper);
  expected_read.manifest_crc64 =
      independent_crc64(manifest_bytes_for_report);
  expected_read.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  expected_read.crc_check_count =
      static_cast<std::uint64_t>(mpi.size()) + 1U;
  expected_read.collective_count = 28U;
  expected_read.rank_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.manifest_crc =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.exact_size_eof =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.fingerprint =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.partition =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.transaction_entry =
      hundun::flow::CheckpointV2CheckStatus::passed;
  expected_read.publication =
      hundun::flow::CheckpointV2CheckStatus::passed;
  require_exact_product_report(
      restored.report(),
      {expected_read, independent_report_fingerprint(expected_read)});
  const auto restored_deep = CheckpointAccess::snapshot(destination);
  HUNDUN_CHECK(!restored_deep.rollback_snapshot_valid);
  HUNDUN_CHECK(!restored_deep.attempt_active);
  HUNDUN_CHECK(!restored_deep.commit_prepared);
  if (!acceptance) {
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(directory);
    return;
  }
  const auto expected_write_failure =
      [&](hundun::flow::CheckpointV2Phase phase, int lowest_failing_rank,
          std::uint64_t file_count, std::uint64_t crc_check_count,
          std::uint64_t collective_count,
          hundun::flow::CheckpointV2CheckStatus rank_crc,
          hundun::flow::CheckpointV2CheckStatus exact_size_eof,
          hundun::flow::CheckpointV2CheckStatus publication) {
        hundun::flow::detail::CheckpointV2ReportValues values;
        values.operation = hundun::flow::CheckpointV2Operation::write;
        values.disposition = hundun::flow::CheckpointV2Disposition::failed;
        values.reason = hundun::flow::CheckpointV2FailureReason::filesystem;
        values.phase = phase;
        values.rank = mpi.rank();
        values.lowest_failing_rank = lowest_failing_rank;
        values.step = metadata.step;
        values.time_s = metadata.time_s;
        values.local_logical_bytes = independent_rank_logical;
        values.local_actual_bytes = expected_rank_actual;
        values.local_crc64 = independent_crc64(actual_rank_wrapper);
        values.file_count = file_count;
        values.crc_check_count = crc_check_count;
        values.collective_count = collective_count;
        values.rank_crc = rank_crc;
        values.exact_size_eof = exact_size_eof;
        values.fingerprint = hundun::flow::CheckpointV2CheckStatus::passed;
        values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
        values.publication = publication;
        return ExpectedProductReport{
            values, independent_report_fingerprint(values)};
      };
  const auto expected_read_file_failure =
      [&](hundun::flow::CheckpointV2Phase phase, int lowest_failing_rank) {
        hundun::flow::detail::CheckpointV2ReportValues values;
        values.operation = hundun::flow::CheckpointV2Operation::read;
        values.disposition = hundun::flow::CheckpointV2Disposition::failed;
        values.reason = hundun::flow::CheckpointV2FailureReason::filesystem;
        values.phase = phase;
        values.rank = mpi.rank();
        values.lowest_failing_rank = lowest_failing_rank;
        values.step = destination_metadata.step;
        values.time_s = destination_metadata.time_s;
        values.transaction_entry =
            hundun::flow::CheckpointV2CheckStatus::passed;
        values.rollback = hundun::flow::CheckpointV2CheckStatus::passed;
        if (phase == hundun::flow::CheckpointV2Phase::marker_read) {
          values.collective_count = 10U;
        } else if (phase ==
                   hundun::flow::CheckpointV2Phase::manifest_read) {
          values.file_count = 1U;
          values.collective_count = 11U;
        } else {
          values.file_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
          values.crc_check_count = static_cast<std::uint64_t>(mpi.size());
          values.collective_count = 26U;
          values.manifest_crc64 =
              independent_crc64(manifest_bytes_for_report);
          values.manifest_crc =
              hundun::flow::CheckpointV2CheckStatus::passed;
          values.fingerprint =
              hundun::flow::CheckpointV2CheckStatus::passed;
          values.partition =
              hundun::flow::CheckpointV2CheckStatus::passed;
          if (mpi.rank() != lowest_failing_rank) {
            values.local_logical_bytes = independent_rank_logical;
            values.local_actual_bytes = expected_rank_actual;
            values.local_crc64 = independent_crc64(actual_rank_wrapper);
            values.rank_crc =
                hundun::flow::CheckpointV2CheckStatus::passed;
          }
        }
        return ExpectedProductReport{
            values, independent_report_fingerprint(values)};
      };
  enum class LateReadEvidence {
    partition,
    global_state,
    physical_state,
    final_success_boundary
  };
  const auto expected_late_read_report =
      [&](const std::filesystem::path &checkpoint_directory,
          const hundun::flow::test::CheckpointV2DeepSnapshot &before,
          LateReadEvidence evidence, int lowest_failing_rank) {
        hundun::flow::detail::CheckpointV2ReportValues values;
        values.operation = hundun::flow::CheckpointV2Operation::read;
        values.disposition = hundun::flow::CheckpointV2Disposition::failed;
        values.reason = evidence == LateReadEvidence::partition
                            ? hundun::flow::CheckpointV2FailureReason::layout
                            : hundun::flow::CheckpointV2FailureReason::state;
        values.phase = evidence == LateReadEvidence::partition
                           ? hundun::flow::CheckpointV2Phase::manifest_read
                           : hundun::flow::CheckpointV2Phase::restore_prepare;
        values.rank = mpi.rank();
        values.lowest_failing_rank = lowest_failing_rank;
        values.step = before.metadata.step;
        values.time_s = before.metadata.time_s;
        values.manifest_crc64 = independent_crc64(
            read_file_bytes(checkpoint_directory / "manifest.v2.bin"));
        values.file_count = 2U;
        values.crc_check_count = 1U;
        values.manifest_crc =
            hundun::flow::CheckpointV2CheckStatus::passed;
        values.exact_size_eof =
            hundun::flow::CheckpointV2CheckStatus::passed;
        values.partition =
            evidence == LateReadEvidence::partition
                ? hundun::flow::CheckpointV2CheckStatus::failed
                : hundun::flow::CheckpointV2CheckStatus::passed;
        values.fingerprint =
            evidence == LateReadEvidence::partition
                ? hundun::flow::CheckpointV2CheckStatus::not_checked
                : hundun::flow::CheckpointV2CheckStatus::passed;
        values.transaction_entry =
            hundun::flow::CheckpointV2CheckStatus::passed;
        values.rollback = hundun::flow::CheckpointV2CheckStatus::passed;
        if (evidence == LateReadEvidence::partition) {
          values.collective_count = 13U;
        } else if (evidence == LateReadEvidence::global_state) {
          values.collective_count = 21U;
        } else {
          const auto rank_bytes =
              read_file_bytes(rank_file_path(checkpoint_directory, mpi.rank()));
          values.local_logical_bytes = independent_rank_logical;
          values.local_actual_bytes = rank_bytes.size();
          values.local_crc64 = independent_crc64(rank_bytes);
          values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
          values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
          values.rank_crc = hundun::flow::CheckpointV2CheckStatus::passed;
          if (evidence == LateReadEvidence::physical_state) {
            values.collective_count = 23U;
          } else {
            values.global_logical_bytes = expected_global_logical;
            values.global_actual_bytes =
                expected_rank_actual_sum + expected_manifest_size + 40U;
            values.collective_count = 29U;
          }
        }
        return ExpectedProductReport{
            values, independent_report_fingerprint(values)};
      };
  if (mpi.size() > 1) {
    using FlowPoint = hundun::flow::test::CheckpointV2PreparationPoint;
    const std::array flow_points{
        FlowPoint::local_layout,      FlowPoint::local_topology,
        FlowPoint::local_geometry,    FlowPoint::topology_common,
        FlowPoint::geometry_common,   FlowPoint::resolved_case,
        FlowPoint::boundary_registry, FlowPoint::field_schema,
        FlowPoint::common_authority};
    for (std::size_t index = 0; index < flow_points.size(); ++index) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-flow-preparation-" + std::to_string(mpi.size()) +
           "-" + std::to_string(index));
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
      mpi.barrier();
      if (mpi.rank() == mpi.size() - 1)
        hundun::flow::test::set_checkpoint_v2_preparation_fault(
            flow_points[index]);
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::state);
      HUNDUN_CHECK(failed.phase() ==
                   hundun::flow::CheckpointV2Phase::preflight);
      HUNDUN_CHECK(failed.lowest_failing_rank() == mpi.size() - 1);
      require_consistent_product_report(mpi, failed);
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory));
    }

    using RuntimePoint = hundun::runtime::checkpoint_v2::test::PreparationPoint;
    const std::array runtime_points{RuntimePoint::opaque_bytes_buffer,
                                    RuntimePoint::allgather_result};
    for (std::size_t index = 0; index < runtime_points.size(); ++index) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-runtime-preparation-" + std::to_string(mpi.size()) +
           "-" + std::to_string(index));
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
      mpi.barrier();
      if (mpi.rank() == mpi.size() - 1)
        hundun::runtime::checkpoint_v2::test::set_preparation_fault(
            runtime_points[index]);
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::state);
      HUNDUN_CHECK(failed.phase() ==
                   hundun::flow::CheckpointV2Phase::preflight);
      HUNDUN_CHECK(failed.lowest_failing_rank() == mpi.size() - 1);
      require_consistent_product_report(mpi, failed);
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory));
    }

    struct LateRuntimePreparation final {
      RuntimePoint point;
      std::uint32_t calls_before;
    };
    const std::array late_runtime_points{
        LateRuntimePreparation{RuntimePoint::allgather_result, 2U},
        LateRuntimePreparation{RuntimePoint::allreduce_workspace, 0U}};
    for (std::size_t index = 0; index < late_runtime_points.size(); ++index) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-late-runtime-preparation-" +
           std::to_string(mpi.size()) + "-" + std::to_string(index));
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
      mpi.barrier();
      if (mpi.rank() == mpi.size() - 1)
        hundun::runtime::checkpoint_v2::test::set_preparation_fault(
            late_runtime_points[index].point,
            late_runtime_points[index].calls_before);
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::state);
      HUNDUN_CHECK(failed.phase() ==
                   hundun::flow::CheckpointV2Phase::rank_publish);
      HUNDUN_CHECK(failed.lowest_failing_rank() == mpi.size() - 1);
      require_consistent_product_report(mpi, failed);
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory / "COMPLETED"));
      mpi.barrier();
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
    }

    {
      auto failed_destination = make_destination_state();
      failed_destination.seed_accepted_layers(different, different);
      const auto before = CheckpointAccess::snapshot(failed_destination);
      if (mpi.rank() == mpi.size() - 1)
        hundun::flow::test::set_checkpoint_v2_preparation_fault(
            FlowPoint::local_topology);
      const auto failed = hundun::flow::read_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config,
          failed_destination, directory);
      HUNDUN_CHECK(!failed.restored());
      HUNDUN_CHECK(failed.report().reason() ==
                   hundun::flow::CheckpointV2FailureReason::state);
      HUNDUN_CHECK(failed.report().phase() ==
                   hundun::flow::CheckpointV2Phase::restore_prepare);
      HUNDUN_CHECK(failed.report().lowest_failing_rank() == mpi.size() - 1);
      require_destination_report_authority(mpi, before, failed.report());
      HUNDUN_CHECK(
          hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
              before, CheckpointAccess::snapshot(failed_destination)));
    }

    {
      auto failed_destination = make_destination_state();
      failed_destination.seed_accepted_layers(different, different);
      const auto before = CheckpointAccess::snapshot(failed_destination);
      if (mpi.rank() == mpi.size() - 1)
        hundun::runtime::checkpoint_v2::test::set_preparation_fault(
            RuntimePoint::allreduce_workspace);
      const auto failed = hundun::flow::read_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config,
          failed_destination, directory);
      HUNDUN_CHECK(!failed.restored());
      HUNDUN_CHECK(failed.report().reason() ==
                   hundun::flow::CheckpointV2FailureReason::state);
      HUNDUN_CHECK(failed.report().phase() ==
                   hundun::flow::CheckpointV2Phase::restore_prepare);
      HUNDUN_CHECK(failed.report().lowest_failing_rank() == mpi.size() - 1);
      require_destination_report_authority(mpi, before, failed.report());
      HUNDUN_CHECK(
          hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
              before, CheckpointAccess::snapshot(failed_destination)));
    }
  }

  for (int failing_rank = 0; failing_rank < mpi.size(); ++failing_rank) {
    auto failed_destination = make_destination_state();
    failed_destination.seed_accepted_layers(different, different);
    const auto before = CheckpointAccess::snapshot(failed_destination);
    const auto old_views = CheckpointAccess::density_views(failed_destination);
    if (mpi.rank() == failing_rank)
      hundun::flow::test::set_checkpoint_v2_preparation_fault(
          hundun::flow::test::CheckpointV2PreparationPoint::
              final_success_boundary);
    const auto failed = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config,
        failed_destination, directory);
    HUNDUN_CHECK(!failed.restored());
    HUNDUN_CHECK(failed.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(failed.report().phase() ==
                 hundun::flow::CheckpointV2Phase::restore_prepare);
    HUNDUN_CHECK(failed.report().lowest_failing_rank() == failing_rank);
    HUNDUN_CHECK(failed.report().publication_status() ==
                 hundun::flow::CheckpointV2CheckStatus::not_checked);
    require_destination_report_authority(mpi, before, failed.report());
    require_exact_product_report(
        failed.report(),
        expected_late_read_report(
            directory, before, LateReadEvidence::final_success_boundary,
            failing_rank));
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        before, CheckpointAccess::snapshot(failed_destination)));
    for (const auto &view : old_views)
      HUNDUN_CHECK(rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));
  }

  using NumericPoint = hundun::runtime::checkpoint_v2::test::NumericFilePoint;
  const std::array rank_write_points{
      NumericPoint::write_status, NumericPoint::write_open,
      NumericPoint::write_body,   NumericPoint::write_flush,
      NumericPoint::write_close,  NumericPoint::read_status,
      NumericPoint::read_size,    NumericPoint::read_open,
      NumericPoint::read_body,    NumericPoint::read_close};
  for (int failing_rank = 0; failing_rank < mpi.size(); ++failing_rank) {
    for (std::size_t point_index = 0; point_index < rank_write_points.size();
         ++point_index) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-rank-file-phase-" + std::to_string(mpi.size()) + "-" +
           std::to_string(failing_rank) + "-" + std::to_string(point_index));
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
      mpi.barrier();
      if (mpi.rank() == failing_rank)
        hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(
            rank_write_points[point_index]);
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::filesystem);
      HUNDUN_CHECK(failed.phase() ==
                   hundun::flow::CheckpointV2Phase::rank_temporary_file);
      HUNDUN_CHECK(failed.lowest_failing_rank() == failing_rank);
      require_consistent_product_report(mpi, failed);
      require_exact_product_report(
          failed,
          expected_write_failure(
              hundun::flow::CheckpointV2Phase::rank_temporary_file,
              failing_rank, 0U,
              static_cast<std::uint64_t>(mpi.size() - 1), 22U,
              mpi.rank() == failing_rank
                  ? hundun::flow::CheckpointV2CheckStatus::not_checked
                  : hundun::flow::CheckpointV2CheckStatus::passed,
              hundun::flow::CheckpointV2CheckStatus::not_checked,
              hundun::flow::CheckpointV2CheckStatus::not_checked));
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory / "COMPLETED"));
      mpi.barrier();
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
    }
    for (const auto point :
         {NumericPoint::publish_status, NumericPoint::publish_rename}) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-rank-publish-phase-" + std::to_string(mpi.size()) +
           "-" + std::to_string(failing_rank) + "-" +
           std::to_string(static_cast<unsigned>(point)));
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
      mpi.barrier();
      if (mpi.rank() == failing_rank)
        hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(point);
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::filesystem);
      HUNDUN_CHECK(failed.phase() ==
                   hundun::flow::CheckpointV2Phase::rank_publish);
      HUNDUN_CHECK(failed.lowest_failing_rank() == failing_rank);
      require_consistent_product_report(mpi, failed);
      require_exact_product_report(
          failed,
          expected_write_failure(
              hundun::flow::CheckpointV2Phase::rank_publish, failing_rank,
              static_cast<std::uint64_t>(mpi.size() - 1),
              static_cast<std::uint64_t>(mpi.size()), 23U,
              hundun::flow::CheckpointV2CheckStatus::passed,
              hundun::flow::CheckpointV2CheckStatus::passed,
              hundun::flow::CheckpointV2CheckStatus::failed));
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory / "COMPLETED"));
      mpi.barrier();
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
    }
  }

  for (const auto point :
       {NumericPoint::directory_status, NumericPoint::parent_status,
        NumericPoint::directory_create}) {
    const auto failure_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-directory-phase-" + std::to_string(mpi.size()) + "-" +
         std::to_string(static_cast<unsigned>(point)));
    if (mpi.rank() == 0) {
      std::filesystem::remove_all(failure_directory);
      hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(point);
    }
    mpi.barrier();
    const auto failed = hundun::flow::write_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config, source,
        controller_state, std::nullopt, failure_directory);
    HUNDUN_CHECK(failed.disposition() ==
                 hundun::flow::CheckpointV2Disposition::failed);
    HUNDUN_CHECK(failed.reason() ==
                 hundun::flow::CheckpointV2FailureReason::filesystem);
    HUNDUN_CHECK(failed.phase() ==
                 hundun::flow::CheckpointV2Phase::rank_temporary_file);
    HUNDUN_CHECK(failed.lowest_failing_rank() == 0);
    require_consistent_product_report(mpi, failed);
    require_exact_product_report(
        failed,
        expected_write_failure(
            hundun::flow::CheckpointV2Phase::rank_temporary_file, 0, 0U, 0U,
            20U, hundun::flow::CheckpointV2CheckStatus::not_checked,
            hundun::flow::CheckpointV2CheckStatus::not_checked,
            hundun::flow::CheckpointV2CheckStatus::not_checked));
    HUNDUN_CHECK(!std::filesystem::exists(failure_directory));
  }

  for (const auto &phase_case :
       {std::pair{1U, hundun::flow::CheckpointV2Phase::manifest},
        std::pair{2U, hundun::flow::CheckpointV2Phase::completed_marker}}) {
    for (const auto point : rank_write_points) {
      const auto failure_directory =
          std::filesystem::temp_directory_path() /
          ("hundun-task23-authority-write-phase-" + std::to_string(mpi.size()) +
           "-" + std::to_string(phase_case.first) + "-" +
           std::to_string(static_cast<unsigned>(point)));
      if (mpi.rank() == 0) {
        std::filesystem::remove_all(failure_directory);
        hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(
            point, phase_case.first);
      }
      mpi.barrier();
      const auto failed = hundun::flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, source,
          controller_state, std::nullopt, failure_directory);
      HUNDUN_CHECK(failed.disposition() ==
                   hundun::flow::CheckpointV2Disposition::failed);
      HUNDUN_CHECK(failed.reason() ==
                   hundun::flow::CheckpointV2FailureReason::filesystem);
      HUNDUN_CHECK(failed.phase() == phase_case.second);
      HUNDUN_CHECK(failed.lowest_failing_rank() == 0);
      require_consistent_product_report(mpi, failed);
      auto expected_late_failure =
          expected_write_failure(
              phase_case.second, 0,
              phase_case.second == hundun::flow::CheckpointV2Phase::manifest
                  ? static_cast<std::uint64_t>(mpi.size())
                  : static_cast<std::uint64_t>(mpi.size()) + 1U,
              phase_case.second == hundun::flow::CheckpointV2Phase::manifest
                  ? static_cast<std::uint64_t>(mpi.size())
                  : static_cast<std::uint64_t>(mpi.size()) + 1U,
              phase_case.second == hundun::flow::CheckpointV2Phase::manifest
                  ? 30U
                  : 33U,
              hundun::flow::CheckpointV2CheckStatus::passed,
              hundun::flow::CheckpointV2CheckStatus::passed,
              hundun::flow::CheckpointV2CheckStatus::failed);
      if (phase_case.second ==
          hundun::flow::CheckpointV2Phase::completed_marker) {
        expected_late_failure.values.global_logical_bytes =
            expected_global_logical;
        expected_late_failure.values.global_actual_bytes =
            expected_rank_actual_sum + expected_manifest_size + 40U;
        expected_late_failure.values.manifest_crc64 =
            independent_crc64(manifest_bytes_for_report);
        expected_late_failure.values.manifest_crc =
            hundun::flow::CheckpointV2CheckStatus::passed;
      }
      expected_late_failure.semantic_fingerprint =
          independent_report_fingerprint(expected_late_failure.values);
      require_exact_product_report(failed, expected_late_failure);
      HUNDUN_CHECK(!std::filesystem::exists(failure_directory / "COMPLETED"));
      mpi.barrier();
      if (mpi.rank() == 0)
        std::filesystem::remove_all(failure_directory);
    }
  }

  using ConfigMutation = std::function<void(hundun::config::FlowCaseConfig &)>;
  std::vector<ConfigMutation> resolved_tuple_mutations;
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.schema_version += 1; });
  resolved_tuple_mutations.push_back([](auto &item) {
    item.simulation_type = hundun::config::SimulationType::passive_scalar;
  });
  resolved_tuple_mutations.push_back([](auto &item) {
    item.density_model = hundun::config::DensityModel::constant;
  });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.resources.expected_ranks = 17; });
  resolved_tuple_mutations.push_back([rank_count = mpi.size()](auto &item) {
    item.resources.process_grid = hundun::runtime::Int3{1, 1, rank_count + 1};
  });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.mesh.cells.x += 1; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.mesh.origin_m.x += 0.125; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.mesh.length_m.y += 0.125; });
  resolved_tuple_mutations.push_back([](auto &item) {
    item.mesh.mapping = hundun::config::MeshMapping::analytic_warped_box;
  });
  resolved_tuple_mutations.push_back([](auto &item) {
    item.mesh.warp_amplitude = {0.01, 0.0, 0.0};
  });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.mode = hundun::config::TimeMode::adaptive; });
  resolved_tuple_mutations.push_back([](auto &item) { item.time.steps += 1; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.initial_dt_s *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.min_dt_s *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.max_dt_s *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.cfl_target *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.diffusion_number_target *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.growth_factor *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.retry_factor *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.time.max_retries += 1; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.rho_ref_kg_per_m3 *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.dynamic_viscosity_pa_s *= 1.01; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.inlet_consistency_rtol *= 2.0; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.cp_J_per_kg_K = 1000.0; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.gas_constant_J_per_kg_K = 287.05; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.physics.thermodynamic_pressure_pa = 101325.0; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.scalars[0].name += "-changed"; });
  resolved_tuple_mutations.push_back(
      [](auto &item) { item.scalars[0].diffusivity_m2_per_s = 0.25; });
  for (std::size_t patch = 0; patch < config.boundaries.size(); ++patch) {
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].patch =
          static_cast<hundun::config::PatchName>((patch + 1U) % 6U);
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].type = hundun::config::BoundaryType::no_slip_wall;
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].velocity_m_per_s =
          hundun::runtime::Real3{1.0, 2.0, 3.0};
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].thermal_authority =
          hundun::config::InletThermalAuthority::temperature;
    });
    resolved_tuple_mutations.push_back(
        [patch](auto &item) { item.boundaries[patch].temperature_K = 300.0; });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].enthalpy_J_per_kg = 300000.0;
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].density_kg_per_m3 = 1.2;
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].scalar_values =
          std::vector<hundun::config::InletScalarValue>{{"alpha", 0.2}};
    });
    resolved_tuple_mutations.push_back([patch](auto &item) {
      item.boundaries[patch].pressure_perturbation_pa = 2.0;
    });
  }
  const auto baseline_resolved_fingerprint =
      independent_resolved_fingerprint(config);
  const auto is_constructible_resolved_fingerprint_case =
      [](std::size_t index) {
        return index == 2U || index == 6U || index == 7U ||
               (index >= 10U && index <= 22U) || index == 26U ||
               index == 27U;
      };
  for (std::size_t mutation_index = 0U;
       mutation_index < resolved_tuple_mutations.size(); ++mutation_index) {
    auto changed = config;
    resolved_tuple_mutations[mutation_index](changed);
    HUNDUN_CHECK(independent_resolved_fingerprint(changed) !=
                 baseline_resolved_fingerprint);
    auto fingerprint_destination = make_destination_state();
    fingerprint_destination.seed_accepted_layers(different, different);
    const auto before = CheckpointAccess::snapshot(fingerprint_destination);
    const auto rejected = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, changed,
        fingerprint_destination, directory);
    if (rejected.restored())
      throw std::runtime_error(
          "Task23 resolved tuple mutation was not rejected: " +
          std::to_string(mutation_index));
    if (is_constructible_resolved_fingerprint_case(mutation_index)) {
      require_constructible_fingerprint_failure(
          mpi, directory, before, fingerprint_destination, rejected);
      continue;
    }
    const auto after = CheckpointAccess::snapshot(fingerprint_destination);
    const bool rejected_before_transaction =
        rejected.report().phase() ==
            hundun::flow::CheckpointV2Phase::preflight ||
        rejected.report().phase() ==
            hundun::flow::CheckpointV2Phase::transaction_entry;
    if (rejected_before_transaction) {
      require_destination_report_authority(mpi, before, rejected.report());
      if (mutation_index == 0U)
        require_exact_product_report(
            rejected.report(),
            expected_read_preflight_failure_report(mpi, before));
      HUNDUN_CHECK(
          hundun::flow::test::checkpoint_v2_deep_snapshot_equal(before, after));
    } else {
      HUNDUN_CHECK(
          hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
              before, after));
    }
  }
  const auto global_crc = independent_crc64(independent_global);
  require_authenticated_mutation_matrix(
      independent_global, [&](const Bytes &candidate) {
        return hundun::flow::test::
            checkpoint_v2_authenticate_global_payload_for_test(
                candidate, global_crc, independent_global.size());
      });
  for (const auto offset : {std::size_t{0U}, std::size_t{36U}}) {
    auto rebuilt_crc_invalid_global = independent_global;
    rebuilt_crc_invalid_global[offset] =
        offset == 36U ? 2U : rebuilt_crc_invalid_global[offset] ^ 1U;
    HUNDUN_CHECK(!hundun::flow::test::
          checkpoint_v2_authenticate_global_payload_for_test(
              rebuilt_crc_invalid_global,
              independent_crc64(rebuilt_crc_invalid_global),
              independent_global.size()));
  }
  auto rebuilt_crc_trailing_global = independent_global;
  rebuilt_crc_trailing_global.push_back(0U);
  HUNDUN_CHECK(!hundun::flow::test::
        checkpoint_v2_authenticate_global_payload_for_test(
            rebuilt_crc_trailing_global,
            independent_crc64(rebuilt_crc_trailing_global),
            independent_global.size()));
  const auto rank_payload_crc = independent_crc64(independent_rank);
  require_authenticated_mutation_matrix(
      independent_rank, [&](const Bytes &candidate) {
        return hundun::flow::test::
            checkpoint_v2_authenticate_rank_payload_for_test(
                candidate, source, rank_payload_crc);
      });
  for (const auto offset :
       {std::size_t{0U}, std::size_t{4U}, std::size_t{16U}}) {
    auto rebuilt_crc_invalid_rank = independent_rank;
    rebuilt_crc_invalid_rank[offset] ^= 1U;
    HUNDUN_CHECK(!hundun::flow::test::
          checkpoint_v2_authenticate_rank_payload_for_test(
              rebuilt_crc_invalid_rank, source,
              independent_crc64(rebuilt_crc_invalid_rank)));
  }
  auto rebuilt_crc_trailing_rank = independent_rank;
  rebuilt_crc_trailing_rank.push_back(0U);
  HUNDUN_CHECK(!hundun::flow::test::
        checkpoint_v2_authenticate_rank_payload_for_test(
            rebuilt_crc_trailing_rank, source,
            independent_crc64(rebuilt_crc_trailing_rank)));

  HUNDUN_CHECK(hundun::flow::test::
          checkpoint_v2_deep_snapshot_oracle_is_mutation_sensitive(
              destination));
  for (std::size_t generation_layer = 0; generation_layer < 4U;
       ++generation_layer) {
    auto generation_state = make_state();
    generation_state.seed_accepted_layers(different, different);
    CheckpointAccess::force_generation(
        generation_state, generation_layer,
        std::numeric_limits<std::uint64_t>::max());
    const auto generation_before = CheckpointAccess::snapshot(generation_state);
    const auto generation_views =
        CheckpointAccess::density_views(generation_state);
    const auto rejected_entry = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config,
        generation_state, directory);
    HUNDUN_CHECK(!rejected_entry.restored());
    HUNDUN_CHECK(rejected_entry.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(rejected_entry.report().phase() ==
                 hundun::flow::CheckpointV2Phase::transaction_entry);
    HUNDUN_CHECK(rejected_entry.report().transaction_entry_status() ==
        hundun::flow::CheckpointV2CheckStatus::failed);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_deep_snapshot_equal(
        generation_before, CheckpointAccess::snapshot(generation_state)));
    for (const auto &view : generation_views)
      HUNDUN_CHECK(std::isfinite(view(0, 0, 0, 0)));
  }

  if (mpi.size() > 1) {
    auto partition_config = config;
    partition_config.resources.process_grid =
        hundun::runtime::Int3{1, mpi.size(), 1};
    auto partition_decomposition =
        hundun::runtime::StructuredDecomposition::create(
            mpi, partition_config.mesh.cells, {true, true, true},
            hundun::runtime::DecompositionOptions{
                *partition_config.resources.process_grid});
    hundun::mesh::MeshTopology partition_topology(partition_decomposition);
    hundun::mesh::MeshGeometry partition_geometry(
        partition_topology,
        hundun::mesh::UniformBoxMapping(partition_config.mesh.origin_m,
                                        partition_config.mesh.length_m));
    auto partition_boundaries = hundun::boundary::BoundaryRegistry::create(
        partition_config, partition_topology);
    auto partition_state =
        hundun::flow::FlowState::create(registry,
        {partition_decomposition.local_extent(),
         partition_topology.local_face_count()},
        fields, destination_metadata);
    const auto partition_values = [&] {
          auto values = history;
          const auto partition_cells = partition_topology.owned_cell_count();
          values.density.assign(partition_cells, 2.0);
          values.velocity.assign(partition_cells * 3U, 0.0);
          values.mechanical_pressure.assign(partition_cells, 0.0);
      values.face_velocity.assign(partition_topology.local_face_count() * 3U,
                                  0.0);
      values.face_mass_flux.assign(partition_topology.local_face_count(), 0.0);
          values.transported_cell_fields = {
          std::vector<double>(partition_cells, 600000.0),
              std::vector<double>(partition_cells, 0.4)};
          return values;
        }();
    partition_state.seed_accepted_layers(partition_values, partition_values);
    const auto partition_before = CheckpointAccess::snapshot(partition_state);
    const auto partition_views =
        CheckpointAccess::density_views(partition_state);
    const auto partition_read = hundun::flow::read_checkpoint_v2(
        mpi, partition_decomposition, partition_topology, partition_geometry,
        partition_boundaries, partition_config, partition_state, directory);
    HUNDUN_CHECK(!partition_read.restored());
    HUNDUN_CHECK(partition_read.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::layout);
    HUNDUN_CHECK(partition_read.report().partition_status() ==
                 hundun::flow::CheckpointV2CheckStatus::failed);
    require_destination_report_authority(mpi, partition_before,
                                         partition_read.report());
    require_exact_product_report(
        partition_read.report(),
        expected_late_read_report(directory, partition_before,
                                  LateReadEvidence::partition, 0));
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
            partition_before, CheckpointAccess::snapshot(partition_state)));
    for (const auto &view : partition_views)
      HUNDUN_CHECK(rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));
  }

  hundun::mesh::MeshGeometry changed_geometry(
      topology, hundun::mesh::AnalyticWarpedBoxMapping(config.mesh.origin_m,
                                                       config.mesh.length_m,
                                                       {0.01, -0.005, 0.0025}));
  auto geometry_destination = make_destination_state();
  geometry_destination.seed_accepted_layers(different, different);
  const auto geometry_before = CheckpointAccess::snapshot(geometry_destination);
  const auto geometry_rejected = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, changed_geometry, boundaries, config,
      geometry_destination, directory);
  HUNDUN_CHECK(!geometry_rejected.restored());
  HUNDUN_CHECK(geometry_rejected.report().reason() ==
               hundun::flow::CheckpointV2FailureReason::file_integrity);
  HUNDUN_CHECK(geometry_rejected.report().fingerprint_status() ==
               hundun::flow::CheckpointV2CheckStatus::failed);
  require_destination_report_authority(mpi, geometry_before,
                                       geometry_rejected.report());
  HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
      geometry_before, CheckpointAccess::snapshot(geometry_destination)));

  const auto source_view = source.layer(hundun::flow::FlowLayer::committed)
          .view<double>(fields.density);
  const double source_first = source_view(0, 0, 0, 0);
  const auto duplicate_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller_state, std::nullopt, directory);
  HUNDUN_CHECK(duplicate_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::failed);
  HUNDUN_CHECK(source_view(0, 0, 0, 0) == source_first);

  const std::array read_points{NumericPoint::read_status,
                               NumericPoint::read_size, NumericPoint::read_open,
                               NumericPoint::read_body,
                               NumericPoint::read_close};
  struct NumericReadPhaseCase final {
    std::uint32_t calls_before;
    hundun::flow::CheckpointV2Phase phase;
  };
  const std::array numeric_read_phase_cases{
      NumericReadPhaseCase{0U, hundun::flow::CheckpointV2Phase::marker_read},
      NumericReadPhaseCase{1U, hundun::flow::CheckpointV2Phase::manifest_read},
      NumericReadPhaseCase{2U, hundun::flow::CheckpointV2Phase::rank_read}};
  for (int failing_rank = 0; failing_rank < mpi.size(); ++failing_rank) {
    for (const auto phase_case : numeric_read_phase_cases) {
      for (const auto point : read_points) {
        auto phase_destination = make_destination_state();
        phase_destination.seed_accepted_layers(different, different);
        const auto before = CheckpointAccess::snapshot(phase_destination);
        if (mpi.rank() == failing_rank)
          hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(
              point, phase_case.calls_before);
        const auto failed = hundun::flow::read_checkpoint_v2(
            mpi, decomposition, topology, geometry, boundaries, config,
            phase_destination, directory);
        HUNDUN_CHECK(!failed.restored());
        HUNDUN_CHECK(failed.report().reason() ==
                     hundun::flow::CheckpointV2FailureReason::filesystem);
        HUNDUN_CHECK(failed.report().phase() == phase_case.phase);
        HUNDUN_CHECK(failed.report().lowest_failing_rank() == failing_rank);
        require_destination_report_authority(mpi, before, failed.report());
        require_exact_product_report(
            failed.report(),
            expected_read_file_failure(phase_case.phase, failing_rank));
        HUNDUN_CHECK(
            hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
                before, CheckpointAccess::snapshot(phase_destination)));
      }
    }
    auto inventory_destination = make_destination_state();
    inventory_destination.seed_accepted_layers(different, different);
    const auto inventory_before =
        CheckpointAccess::snapshot(inventory_destination);
    if (mpi.rank() == failing_rank)
      hundun::runtime::checkpoint_v2::test::set_numeric_file_fault(
          NumericPoint::inventory_iteration);
    const auto inventory_failed = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config,
        inventory_destination, directory);
    HUNDUN_CHECK(!inventory_failed.restored());
    HUNDUN_CHECK(inventory_failed.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::filesystem);
    HUNDUN_CHECK(inventory_failed.report().phase() ==
                 hundun::flow::CheckpointV2Phase::manifest_read);
    HUNDUN_CHECK(inventory_failed.report().lowest_failing_rank() ==
                 failing_rank);
    require_destination_report_authority(mpi, inventory_before,
                                         inventory_failed.report());
    auto expected_inventory = expected_read_file_failure(
        hundun::flow::CheckpointV2Phase::manifest_read, failing_rank);
    expected_inventory.values.file_count = 2U;
    expected_inventory.values.crc_check_count = 1U;
    expected_inventory.values.collective_count = 12U;
    expected_inventory.values.manifest_crc64 =
        independent_crc64(manifest_bytes_for_report);
    expected_inventory.values.manifest_crc =
        hundun::flow::CheckpointV2CheckStatus::passed;
    expected_inventory.values.exact_size_eof =
        hundun::flow::CheckpointV2CheckStatus::passed;
    expected_inventory.semantic_fingerprint =
        independent_report_fingerprint(expected_inventory.values);
    require_exact_product_report(inventory_failed.report(),
                                 expected_inventory);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        inventory_before, CheckpointAccess::snapshot(inventory_destination)));
  }

  CheckpointAccess::set_committed_density_ghost(destination, -0.0);
  CheckpointAccess::set_rollback_density_ghost(destination, 7.25);
  const auto before_deep = CheckpointAccess::snapshot(destination);
  const auto old_views = CheckpointAccess::density_views(destination);
  mpi.barrier();
  if (mpi.rank() == 0) {
    const auto manifest = directory / "manifest.v2.bin";
    std::fstream stream(manifest,
                        std::ios::in | std::ios::out | std::ios::binary);
    HUNDUN_CHECK(static_cast<bool>(stream));
    stream.seekg(24);
    char byte{};
    stream.read(&byte, 1);
    HUNDUN_CHECK(static_cast<bool>(stream));
    byte ^= 1;
    stream.seekp(24);
    stream.write(&byte, 1);
    stream.flush();
    HUNDUN_CHECK(static_cast<bool>(stream));
  }
  mpi.barrier();
  const auto failed_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, destination,
      directory);
  HUNDUN_CHECK(!failed_read.restored());
  HUNDUN_CHECK(failed_read.report().reason() ==
               hundun::flow::CheckpointV2FailureReason::file_integrity);
  HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
          before_deep, CheckpointAccess::snapshot(destination)));
  for (const auto &view : old_views)
    HUNDUN_CHECK(rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));

  auto ideal_config = config;
  ideal_config.density_model = hundun::config::DensityModel::ideal_gas;
  ideal_config.scalars = {{"alpha", 0.0}};
  ideal_config.physics.cp_J_per_kg_K = 1000.0;
  ideal_config.physics.gas_constant_J_per_kg_K = 287.05;
  ideal_config.physics.thermodynamic_pressure_pa = 101325.0;
  auto ideal_boundaries =
      hundun::boundary::BoundaryRegistry::create(ideal_config, topology);
  hundun::runtime::FieldRegistry ideal_registry;
  hundun::flow::FlowFieldIds ideal_fields;
  ideal_fields.density =
      ideal_registry.declare_field(physical_cell("rho", "kg/m3", true));
  ideal_fields.velocity = ideal_registry.declare_field(cell("u", 3U, false));
  ideal_fields.mechanical_pressure =
      ideal_registry.declare_field(cell("pi", 1U, false));
  ideal_fields.face_velocity = ideal_registry.declare_field(face("uf", 3U));
  ideal_fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(ideal_registry);
  const auto rho_h =
      ideal_registry.declare_field(physical_cell("rho_h", "J/m3", true));
  const auto ideal_rho_alpha =
      ideal_registry.declare_field(physical_cell("rho_alpha", "kg/m3", true));
  ideal_fields.transported_cell_fields = {rho_h, ideal_rho_alpha};
  ideal_registry.freeze();
  auto ideal_state = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()}, ideal_fields,
      metadata);
  constexpr double temperature = 300.0;
  const double ideal_density = 101325.0 / (287.05 * temperature);
  hundun::flow::FlowLayerValues ideal_values;
  ideal_values.density.resize(cells);
  for (std::size_t cell_id = 0; cell_id < cells; ++cell_id) {
    const double heated_temperature =
        temperature +
        20.0 * static_cast<double>(topology.global_cell_id(cell_id) % 5U);
    ideal_values.density[cell_id] = 101325.0 / (287.05 * heated_temperature);
  }
  ideal_values.velocity.assign(cells * 3U, 0.0);
  ideal_values.mechanical_pressure.assign(cells, 0.0);
  ideal_values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  ideal_values.face_mass_flux.assign(topology.local_face_count(), 0.0);
  ideal_values.transported_cell_fields = {
      std::vector<double>(cells, 101325.0 * 1000.0 / 287.05),
      std::vector<double>(cells)};
  for (std::size_t cell_id = 0; cell_id < cells; ++cell_id)
    ideal_values.transported_cell_fields[1][cell_id] =
        ideal_values.density[cell_id] * 0.2;
  ideal_state.seed_accepted_layers(ideal_values, ideal_values);
  const hundun::flow::IdealGasClosureSpec ideal_spec{rho_h, 1000.0, 287.05,
                                                     101325.0};
  auto initial_closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_state, ideal_spec);
  auto persisted_closure = initial_closure.state();
  persisted_closure.revision = 0U;
  if (mpi.size() > 1) {
    for (int failure_rank = 0; failure_rank < mpi.size(); ++failure_rank) {
      const auto restore_before = CheckpointAccess::snapshot(ideal_state);
      hundun::flow::test::set_ideal_gas_restore_snapshot_preparation_fault(
          failure_rank);
      bool exact_failure = false;
      try {
        static_cast<void>(hundun::flow::IdealGasClosure::restore(
            topology, geometry, ideal_boundaries, mpi, ideal_registry,
            ideal_fields, ideal_state, ideal_spec, persisted_closure));
      } catch (const hundun::runtime::Error &error) {
        exact_failure =
            std::string_view(error.what())
                    .find("ideal-gas closure restore snapshot preparation "
                          "failed") != std::string_view::npos &&
            std::string_view(error.what()).find(std::to_string(failure_rank)) !=
                std::string_view::npos;
      }
      HUNDUN_CHECK(exact_failure);
      HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_deep_snapshot_equal(
          restore_before, CheckpointAccess::snapshot(ideal_state)));
    }
  }
  auto restored_closure = hundun::flow::IdealGasClosure::restore(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_state, ideal_spec, persisted_closure);
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      restored_closure.state(), persisted_closure));
  auto invalid_closure = persisted_closure;
  invalid_closure.revision = std::numeric_limits<std::uint64_t>::max();
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(hundun::flow::IdealGasClosure::restore(
        topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
        ideal_state, ideal_spec, invalid_closure));
  }));

  auto ideal_controller = hundun::flow::Bdf2RetryController::create(
      ideal_config.time, ideal_config.density_model, topology, geometry, mpi,
      ideal_state);
  const auto ideal_controller_state = ideal_controller.state();
  const auto ideal_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-ideal-checkpoint-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(ideal_directory);
  }
  mpi.barrier();
  const auto ideal_written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      ideal_state, ideal_controller_state, persisted_closure, ideal_directory);
  HUNDUN_CHECK(ideal_written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  const auto ideal_max_revision_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-ideal-max-revision-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(ideal_max_revision_directory);
    std::filesystem::copy(ideal_directory, ideal_max_revision_directory,
                          std::filesystem::copy_options::recursive);
  }
  mpi.barrier();
  rebuild_checkpoint_global_u64(mpi, ideal_max_revision_directory, 124U,
                                std::numeric_limits<std::uint64_t>::max());
  auto max_revision_destination = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()}, ideal_fields,
      destination_metadata);
  max_revision_destination.seed_accepted_layers(ideal_values, ideal_values);
  const auto max_revision_before =
      CheckpointAccess::snapshot(max_revision_destination);
  const auto max_revision_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      max_revision_destination, ideal_max_revision_directory);
  HUNDUN_CHECK(!max_revision_read.restored());
  HUNDUN_CHECK(max_revision_read.report().reason() ==
               hundun::flow::CheckpointV2FailureReason::state);
  HUNDUN_CHECK(max_revision_read.report().phase() ==
               hundun::flow::CheckpointV2Phase::restore_prepare);
  require_destination_report_authority(mpi, max_revision_before,
                                       max_revision_read.report());
  HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
      max_revision_before,
      CheckpointAccess::snapshot(max_revision_destination)));
  struct ClosureStateMutation final {
    const char *name;
    std::size_t offset;
    Bytes replacement;
    hundun::flow::CheckpointV2FailureReason reason;
    hundun::flow::CheckpointV2Phase phase;
  };
  const auto closure_u8 = [](std::uint8_t value) {
    Bytes bytes;
    append_u8(bytes, value);
    return bytes;
  };
  const auto closure_f64 = [](double value) {
    Bytes bytes;
    append_f64(bytes, value);
    return bytes;
  };
  const std::array closure_state_mutations{
      ClosureStateMutation{
          "presence", 105U, closure_u8(2U),
          hundun::flow::CheckpointV2FailureReason::file_integrity,
          hundun::flow::CheckpointV2Phase::manifest_read},
      ClosureStateMutation{"mode", 106U, closure_u8(1U),
                           hundun::flow::CheckpointV2FailureReason::state,
                           hundun::flow::CheckpointV2Phase::restore_prepare},
      ClosureStateMutation{"pressure-negative", 107U, closure_f64(-1.0),
                           hundun::flow::CheckpointV2FailureReason::state,
                           hundun::flow::CheckpointV2Phase::restore_prepare},
      ClosureStateMutation{
          "pressure-nan", 107U,
          closure_f64(std::numeric_limits<double>::quiet_NaN()),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      ClosureStateMutation{
          "target-presence", 115U, closure_u8(2U),
          hundun::flow::CheckpointV2FailureReason::file_integrity,
          hundun::flow::CheckpointV2Phase::manifest_read},
      ClosureStateMutation{"target-negative", 116U, closure_f64(-1.0),
                           hundun::flow::CheckpointV2FailureReason::state,
                           hundun::flow::CheckpointV2Phase::restore_prepare},
      ClosureStateMutation{
          "target-nan", 116U,
          closure_f64(std::numeric_limits<double>::quiet_NaN()),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare}};
  for (const auto &mutation : closure_state_mutations) {
    const auto invalid_closure_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-closure-state-" + std::string(mutation.name) + "-" +
         std::to_string(mpi.size()));
    if (mpi.rank() == 0) {
      std::filesystem::remove_all(invalid_closure_directory);
      std::filesystem::copy(ideal_directory, invalid_closure_directory,
                            std::filesystem::copy_options::recursive);
    }
    mpi.barrier();
    rebuild_checkpoint_global_bytes(mpi, invalid_closure_directory,
                                    mutation.offset, mutation.replacement);
    auto invalid_closure_destination = hundun::flow::FlowState::create(
        ideal_registry,
        {decomposition.local_extent(), topology.local_face_count()},
        ideal_fields, metadata);
    invalid_closure_destination.seed_accepted_layers(ideal_values,
                                                     ideal_values);
    const auto invalid_closure_before =
        CheckpointAccess::snapshot(invalid_closure_destination);
    const auto invalid_closure_read = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
        invalid_closure_destination, invalid_closure_directory);
    HUNDUN_CHECK(!invalid_closure_read.restored());
    HUNDUN_CHECK(invalid_closure_read.report().reason() == mutation.reason);
    HUNDUN_CHECK(invalid_closure_read.report().phase() == mutation.phase);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        invalid_closure_before,
        CheckpointAccess::snapshot(invalid_closure_destination)));
    require_consistent_product_report(mpi, invalid_closure_read.report());
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(invalid_closure_directory);
    mpi.barrier();
  }
  auto ideal_destination = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()}, ideal_fields,
      metadata);
  ideal_destination.seed_accepted_layers(ideal_values, ideal_values);
  const auto ideal_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      ideal_destination, ideal_directory);
  HUNDUN_CHECK(ideal_read.restored());
  HUNDUN_CHECK(ideal_read.ideal_gas_closure_state_available());
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      ideal_read.ideal_gas_closure_state(), persisted_closure));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      ideal_read.time_control_state(), ideal_controller_state));
  auto resumed_closure = hundun::flow::IdealGasClosure::restore(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_destination, ideal_spec, ideal_read.ideal_gas_closure_state());
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      resumed_closure.state(), persisted_closure));

  auto ideal_resumed_controller = hundun::flow::Bdf2RetryController::restore(
      ideal_config.time, ideal_config.density_model, topology, geometry, mpi,
      ideal_destination, ideal_read.time_control_state());
  auto ideal_halo = hundun::runtime::HaloExchange::create(
      decomposition, hundun::runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  auto ideal_resumed_halo = hundun::runtime::HaloExchange::create(
      decomposition, hundun::runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  hundun::execution::CpuReferenceContext ideal_execution;
  hundun::execution::CpuReferenceContext ideal_resumed_execution;
  hundun::linear::ConjugateGradientSolver ideal_momentum(ideal_execution, mpi),
      ideal_pressure(ideal_execution, mpi);
  hundun::linear::ConjugateGradientSolver ideal_resumed_momentum(
      ideal_resumed_execution, mpi);
  hundun::linear::ConjugateGradientSolver ideal_resumed_pressure(
      ideal_resumed_execution, mpi);
  hundun::linear::JacobiPreconditioner ideal_mx(ideal_execution),
      ideal_my(ideal_execution), ideal_mz(ideal_execution),
      ideal_pressure_preconditioner(ideal_execution);
  hundun::linear::JacobiPreconditioner ideal_resumed_mx(
      ideal_resumed_execution),
      ideal_resumed_my(ideal_resumed_execution),
      ideal_resumed_mz(ideal_resumed_execution),
      ideal_resumed_pressure_preconditioner(ideal_resumed_execution);
  hundun::flow::MaterialDensityTransportSpec ideal_material_spec;
  ideal_material_spec.enthalpy_density = rho_h;
  ideal_material_spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  ideal_material_spec.scalar_densities = {ideal_rho_alpha};
  ideal_material_spec.scalar_diffusivities_kg_per_m_s = {0.0};
  auto ideal_flow = hundun::flow::FixedStepIdealGasFlow::create(
      decomposition, topology, geometry, ideal_boundaries, mpi, ideal_execution,
      ideal_halo, ideal_momentum, {&ideal_mx, &ideal_my, &ideal_mz},
      ideal_pressure, ideal_pressure_preconditioner, ideal_registry,
      ideal_fields, ideal_material_spec, std::move(restored_closure));
  auto ideal_resumed_flow = hundun::flow::FixedStepIdealGasFlow::create(
      decomposition, topology, geometry, ideal_boundaries, mpi,
      ideal_resumed_execution, ideal_resumed_halo, ideal_resumed_momentum,
      {&ideal_resumed_mx, &ideal_resumed_my, &ideal_resumed_mz},
      ideal_resumed_pressure, ideal_resumed_pressure_preconditioner,
      ideal_registry, ideal_fields, ideal_material_spec,
      std::move(resumed_closure));
  const auto ideal_next =
      ideal_controller.advance(ideal_state, ideal_flow, 0.0, {}, {});
  const auto ideal_resumed_next = ideal_resumed_controller.advance(
      ideal_destination, ideal_resumed_flow, 0.0, {}, {});
  HUNDUN_CHECK(ideal_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(ideal_resumed_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  for (const auto layer :
       {hundun::flow::FlowLayer::history, hundun::flow::FlowLayer::committed,
        hundun::flow::FlowLayer::trial}) {
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        ideal_state.snapshot(layer), ideal_destination.snapshot(layer)));
  }
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      ideal_state.metadata(), ideal_destination.metadata()));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      ideal_controller.state(), ideal_resumed_controller.state()));
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      ideal_flow.closure_state(), ideal_resumed_flow.closure_state()));

  auto open_config = ideal_config;
  open_config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  open_config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
  open_config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::temperature;
  open_config.boundaries[0].temperature_K = 300.0;
  open_config.boundaries[0].enthalpy_J_per_kg = 300000.0;
  open_config.boundaries[0].density_kg_per_m3 = ideal_density;
  open_config.boundaries[0].scalar_values =
      std::vector<hundun::config::InletScalarValue>{{"alpha", 0.2}};
  open_config.boundaries[1].type =
      hundun::config::BoundaryType::pressure_outlet;
  open_config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < open_config.boundaries.size(); ++patch)
    open_config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
  auto open_decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, open_config.mesh.cells, {false, false, false},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology open_topology(open_decomposition);
  hundun::mesh::MeshGeometry open_geometry(
      open_topology, hundun::mesh::UniformBoxMapping(
                         open_config.mesh.origin_m, open_config.mesh.length_m));
  auto open_boundaries =
      hundun::boundary::BoundaryRegistry::create(open_config, open_topology);
  hundun::runtime::FieldRegistry open_registry;
  hundun::flow::FlowFieldIds open_fields;
  open_fields.density =
      open_registry.declare_field(physical_cell("rho", "kg/m3", true));
  open_fields.velocity = open_registry.declare_field(cell("u", 3U, false));
  open_fields.mechanical_pressure =
      open_registry.declare_field(cell("pi", 1U, false));
  open_fields.face_velocity = open_registry.declare_field(face("uf", 3U));
  open_fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(open_registry);
  const auto open_rho_h =
      open_registry.declare_field(physical_cell("rho_h", "J/m3", true));
  const auto open_rho_alpha =
      open_registry.declare_field(physical_cell("rho_alpha", "kg/m3", true));
  open_fields.transported_cell_fields = {open_rho_h, open_rho_alpha};
  open_registry.freeze();
  hundun::flow::FlowLayerValues open_values;
  open_values.density.assign(open_topology.owned_cell_count(), ideal_density);
  open_values.velocity.assign(open_topology.owned_cell_count() * 3U, 0.0);
  for (std::size_t cell_id = 0; cell_id < open_topology.owned_cell_count();
       ++cell_id)
    open_values.velocity[cell_id * 3U] = 1.0;
  open_values.mechanical_pressure.assign(open_topology.owned_cell_count(), 0.0);
  open_values.face_velocity.assign(open_topology.local_face_count() * 3U, 0.0);
  open_values.face_mass_flux.assign(open_topology.local_face_count(), 0.0);
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < open_topology.local_face_count(); ++face_id) {
    const auto area = open_geometry.face_area_vector_m2(
        face_id, hundun::mesh::FaceSide::owner);
    open_values.face_velocity[face_id * 3U] = 1.0;
    open_values.face_mass_flux[face_id] = ideal_density * area.x;
  }
  open_values.transported_cell_fields = {
      std::vector<double>(open_topology.owned_cell_count(),
                          ideal_density * 300000.0),
      std::vector<double>(open_topology.owned_cell_count(),
                          ideal_density * 0.2)};
  const auto make_open_state = [&] {
    auto result = hundun::flow::FlowState::create(
        open_registry,
        {open_decomposition.local_extent(), open_topology.local_face_count()},
        open_fields, metadata);
    result.seed_accepted_layers(open_values, open_values);
    return result;
  };
  const auto make_open_destination_state = [&] {
    auto result = hundun::flow::FlowState::create(
        open_registry,
        {open_decomposition.local_extent(), open_topology.local_face_count()},
        open_fields, destination_metadata);
    result.seed_accepted_layers(open_values, open_values);
    return result;
  };
  auto open_state = make_open_state();
  const hundun::flow::IdealGasClosureSpec open_spec{open_rho_h, 1000.0, 287.05,
                                                    101325.0};
  auto open_closure = hundun::flow::IdealGasClosure::create(
      open_topology, open_geometry, open_boundaries, mpi, open_registry,
      open_fields, open_state, open_spec);
  const auto open_closure_state = open_closure.state();
  HUNDUN_CHECK(open_closure_state.mode ==
               hundun::flow::IdealGasPressureMode::open_fixed);
  HUNDUN_CHECK(!open_closure_state.target_mass_kg.has_value());
  auto open_controller = hundun::flow::Bdf2RetryController::create(
      open_config.time, open_config.density_model, open_topology, open_geometry,
      mpi, open_state);
  const auto open_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-open-ideal-continuation-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(open_directory);
  mpi.barrier();
  const auto open_write = hundun::flow::write_checkpoint_v2(
      mpi, open_decomposition, open_topology, open_geometry, open_boundaries,
      open_config, open_state, open_controller.state(), open_closure_state,
      open_directory);
  HUNDUN_CHECK(open_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  const auto open_boundary_fingerprint =
      independent_boundary_fingerprint(open_boundaries);
  for (std::size_t inlet_patch = 0U;
       inlet_patch < open_config.boundaries.size(); ++inlet_patch) {
    for (std::size_t authority_case = 0U; authority_case < 2U;
         ++authority_case) {
    auto changed_open_config = open_config;
    for (std::size_t patch = 0U;
         patch < changed_open_config.boundaries.size(); ++patch) {
      const auto patch_name = changed_open_config.boundaries[patch].patch;
      changed_open_config.boundaries[patch] = {};
      changed_open_config.boundaries[patch].patch = patch_name;
      changed_open_config.boundaries[patch].type =
          hundun::config::BoundaryType::symmetry;
    }
    const auto outlet_patch = inlet_patch ^ 1U;
    auto &inlet = changed_open_config.boundaries[inlet_patch];
    inlet.type = hundun::config::BoundaryType::velocity_inlet;
      const double speed = 1.0 + 0.01 * static_cast<double>(inlet_patch);
    const std::array<hundun::runtime::Real3, 6> inward_velocity{
        hundun::runtime::Real3{speed, 0.0, 0.0},
        hundun::runtime::Real3{-speed, 0.0, 0.0},
        hundun::runtime::Real3{0.0, speed, 0.0},
        hundun::runtime::Real3{0.0, -speed, 0.0},
        hundun::runtime::Real3{0.0, 0.0, speed},
        hundun::runtime::Real3{0.0, 0.0, -speed}};
    inlet.velocity_m_per_s = inward_velocity[inlet_patch];
    inlet.thermal_authority =
          authority_case == 0U
            ? hundun::config::InletThermalAuthority::temperature
            : hundun::config::InletThermalAuthority::enthalpy;
    inlet.temperature_K = 300.0;
    inlet.enthalpy_J_per_kg = 300000.0;
    inlet.density_kg_per_m3 = ideal_density;
    inlet.scalar_values =
          std::vector<hundun::config::InletScalarValue>{{"alpha", 0.2}};
    auto &outlet = changed_open_config.boundaries[outlet_patch];
    outlet.type = hundun::config::BoundaryType::pressure_outlet;
      outlet.pressure_perturbation_pa = 0.5 + static_cast<double>(inlet_patch);
      auto changed_open_boundaries = hundun::boundary::BoundaryRegistry::create(
          changed_open_config, open_topology);
    HUNDUN_CHECK(independent_boundary_fingerprint(changed_open_boundaries) !=
                 open_boundary_fingerprint);
    auto changed_open_state = make_open_destination_state();
    const auto changed_open_before =
        CheckpointAccess::snapshot(changed_open_state);
      const auto changed_open_views =
          CheckpointAccess::density_views(changed_open_state);
    const auto changed_open_read = hundun::flow::read_checkpoint_v2(
        mpi, open_decomposition, open_topology, open_geometry,
        changed_open_boundaries, changed_open_config, changed_open_state,
        open_directory);
    require_constructible_fingerprint_failure(
        mpi, open_directory, changed_open_before, changed_open_state,
        changed_open_read);
      for (const auto &view : changed_open_views)
        HUNDUN_CHECK(rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));

    const auto variant_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-open-boundary-variant-" +
           std::to_string(inlet_patch) + "-" + std::to_string(authority_case) +
           "-" + std::to_string(mpi.size()));
    if (mpi.rank() == 0)
      std::filesystem::remove_all(variant_directory);
    mpi.barrier();
    const auto variant_write = hundun::flow::write_checkpoint_v2(
        mpi, open_decomposition, open_topology, open_geometry,
        changed_open_boundaries, changed_open_config, open_state,
        open_controller.state(), open_closure_state, variant_directory);
    HUNDUN_CHECK(variant_write.disposition() ==
                 hundun::flow::CheckpointV2Disposition::completed);
    const auto require_single_boundary_mutation =
        [&](const ConfigMutation &mutate) {
          auto single = changed_open_config;
          mutate(single);
            auto single_boundaries = hundun::boundary::BoundaryRegistry::create(
                single, open_topology);
          HUNDUN_CHECK(
              independent_boundary_fingerprint(single_boundaries) !=
                  independent_boundary_fingerprint(changed_open_boundaries) ||
              independent_resolved_fingerprint(single) !=
                  independent_resolved_fingerprint(changed_open_config));
          auto single_state = make_open_destination_state();
          const auto single_before = CheckpointAccess::snapshot(single_state);
            const auto single_views =
                CheckpointAccess::density_views(single_state);
          const auto single_read = hundun::flow::read_checkpoint_v2(
              mpi, open_decomposition, open_topology, open_geometry,
              single_boundaries, single, single_state, variant_directory);
          require_constructible_fingerprint_failure(
              mpi, variant_directory, single_before, single_state,
              single_read);
            for (const auto &view : single_views)
              HUNDUN_CHECK(
                  rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));
        };
    require_single_boundary_mutation([inlet_patch](auto &single) {
        auto &velocity = *single.boundaries[inlet_patch].velocity_m_per_s;
      velocity.x += velocity.x == 0.0 ? 0.0 : 0.001;
      velocity.y += velocity.y == 0.0 ? 0.0 : 0.001;
      velocity.z += velocity.z == 0.0 ? 0.0 : 0.001;
    });
    require_single_boundary_mutation([inlet_patch](auto &single) {
        auto &authority = *single.boundaries[inlet_patch].thermal_authority;
      authority =
          authority == hundun::config::InletThermalAuthority::temperature
              ? hundun::config::InletThermalAuthority::enthalpy
              : hundun::config::InletThermalAuthority::temperature;
    });
    require_single_boundary_mutation([inlet_patch](auto &single) {
      const auto authority =
          *single.boundaries[inlet_patch].thermal_authority;
        if (authority == hundun::config::InletThermalAuthority::temperature)
        single.boundaries[inlet_patch].enthalpy_J_per_kg.reset();
      else
        single.boundaries[inlet_patch].temperature_K.reset();
    });
      const ConfigMutation r5_authoritative_thermal_numeric =
          [inlet_patch](auto &single) {
            constexpr double changed_temperature = 310.0;
            auto &changed = single.boundaries[inlet_patch];
            changed.temperature_K = changed_temperature;
            changed.enthalpy_J_per_kg = 1000.0 * changed_temperature;
            changed.density_kg_per_m3 =
                101325.0 / (287.05 * changed_temperature);
          };
      require_single_boundary_mutation(r5_authoritative_thermal_numeric);
      const ConfigMutation r5_density_optional_presence =
          [inlet_patch](auto &single) {
            single.boundaries[inlet_patch].density_kg_per_m3.reset();
          };
      require_single_boundary_mutation(r5_density_optional_presence);
      const ConfigMutation r5_inlet_scalar_value = [inlet_patch](auto &single) {
        single.boundaries[inlet_patch].scalar_values->front().value += 0.05;
      };
      require_single_boundary_mutation(r5_inlet_scalar_value);
      require_single_boundary_mutation([outlet_patch](auto &single) {
        *single.boundaries[outlet_patch].pressure_perturbation_pa += 0.125;
      });
      for (std::size_t wall_patch = 0U;
           wall_patch < changed_open_config.boundaries.size(); ++wall_patch) {
        if (wall_patch != inlet_patch && wall_patch != outlet_patch)
          require_single_boundary_mutation([wall_patch](auto &single) {
            single.boundaries[wall_patch].type =
                hundun::config::BoundaryType::no_slip_wall;
          });
      }

      {
        const auto require_scalar_schema_mutation =
            [&](hundun::config::FlowCaseConfig single,
                const std::optional<std::string> &field_name) {
              auto single_boundaries =
                  hundun::boundary::BoundaryRegistry::create(single,
                                                             open_topology);
              HUNDUN_CHECK(
                  independent_boundary_fingerprint(single_boundaries) !=
                      independent_boundary_fingerprint(
                          changed_open_boundaries) ||
                  independent_resolved_fingerprint(single) !=
                      independent_resolved_fingerprint(changed_open_config));
              hundun::runtime::FieldRegistry single_registry;
              hundun::flow::FlowFieldIds single_fields;
              single_fields.density = single_registry.declare_field(
                  physical_cell("rho", "kg/m3", true));
              single_fields.velocity =
                  single_registry.declare_field(cell("u", 3U, false));
              single_fields.mechanical_pressure =
                  single_registry.declare_field(cell("pi", 1U, false));
              single_fields.face_velocity =
                  single_registry.declare_field(face("uf", 3U));
              single_fields.face_mass_flux =
                  hundun::finite_volume::declare_face_mass_flux(
                      single_registry);
              single_fields.transported_cell_fields = {
                  single_registry.declare_field(
                      physical_cell("rho_h", "J/m3", true))};
              if (field_name)
                single_fields.transported_cell_fields.push_back(
                    single_registry.declare_field(
                        physical_cell(field_name->c_str(), "kg/m3", true)));
              single_registry.freeze();
              auto single_values = open_values;
              single_values.transported_cell_fields.resize(field_name ? 2U
                                                                      : 1U);
              auto single_state = hundun::flow::FlowState::create(
                  single_registry,
                  {open_decomposition.local_extent(),
                   open_topology.local_face_count()},
                  single_fields, destination_metadata);
              single_state.seed_accepted_layers(single_values, single_values);
              const auto single_before =
                  CheckpointAccess::snapshot(single_state);
              const auto single_views =
                  CheckpointAccess::density_views(single_state);
              const auto single_read = hundun::flow::read_checkpoint_v2(
                  mpi, open_decomposition, open_topology, open_geometry,
                  single_boundaries, single, single_state, variant_directory);
              require_constructible_fingerprint_failure(
                  mpi, variant_directory, single_before, single_state,
                  single_read);
              for (const auto &view : single_views)
                HUNDUN_CHECK(
                    rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));
            };
        auto r5_scalar_list_membership = changed_open_config;
        r5_scalar_list_membership.scalars.clear();
        r5_scalar_list_membership.boundaries[inlet_patch].scalar_values =
            std::vector<hundun::config::InletScalarValue>{};
        require_scalar_schema_mutation(r5_scalar_list_membership, std::nullopt);
        auto r5_scalar_name = changed_open_config;
        r5_scalar_name.scalars.front().name = "beta";
        r5_scalar_name.boundaries[inlet_patch].scalar_values->front().name =
            "beta";
        require_scalar_schema_mutation(r5_scalar_name, std::string("rho_beta"));
      }
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(variant_directory);
  }
  }

  auto material_open_config = open_config;
  material_open_config.density_model = hundun::config::DensityModel::material;
  material_open_config.physics.cp_J_per_kg_K.reset();
  material_open_config.physics.gas_constant_J_per_kg_K.reset();
  material_open_config.physics.thermodynamic_pressure_pa.reset();
  for (auto &boundary : material_open_config.boundaries) {
    const auto patch_name = boundary.patch;
    boundary = {};
    boundary.patch = patch_name;
    boundary.type = hundun::config::BoundaryType::symmetry;
  }
  auto material_open_values = open_values;
  std::fill(material_open_values.density.begin(),
            material_open_values.density.end(), 1.0);
  std::fill(material_open_values.face_mass_flux.begin(),
            material_open_values.face_mass_flux.end(), 0.0);
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < open_topology.local_face_count(); ++face_id) {
    const auto area = open_geometry.face_area_vector_m2(
        face_id, hundun::mesh::FaceSide::owner);
    material_open_values.face_mass_flux[face_id] = area.x;
  }
  std::fill(material_open_values.transported_cell_fields[0].begin(),
            material_open_values.transported_cell_fields[0].end(), 300000.0);
  std::fill(material_open_values.transported_cell_fields[1].begin(),
            material_open_values.transported_cell_fields[1].end(), 0.2);
  const auto make_material_open_state = [&] {
    auto result = hundun::flow::FlowState::create(
        open_registry,
        {open_decomposition.local_extent(), open_topology.local_face_count()},
        open_fields, metadata);
    result.seed_accepted_layers(material_open_values, material_open_values);
    return result;
  };
  const auto make_material_open_destination_state = [&] {
    auto result = hundun::flow::FlowState::create(
        open_registry,
        {open_decomposition.local_extent(), open_topology.local_face_count()},
        open_fields, destination_metadata);
    result.seed_accepted_layers(material_open_values, material_open_values);
    return result;
  };
  const std::array<hundun::runtime::Real3, 6> material_inward_velocity{
      hundun::runtime::Real3{1.0, 0.0, 0.0},
      hundun::runtime::Real3{-1.0, 0.0, 0.0},
      hundun::runtime::Real3{0.0, 1.0, 0.0},
      hundun::runtime::Real3{0.0, -1.0, 0.0},
      hundun::runtime::Real3{0.0, 0.0, 1.0},
      hundun::runtime::Real3{0.0, 0.0, -1.0}};
  for (std::size_t material_inlet_patch = 0U;
       material_inlet_patch < material_open_config.boundaries.size();
       ++material_inlet_patch) {
    auto material_patch_config = material_open_config;
    const auto material_outlet_patch = material_inlet_patch ^ 1U;
    auto &material_inlet =
        material_patch_config.boundaries[material_inlet_patch];
    material_inlet.type = hundun::config::BoundaryType::velocity_inlet;
    material_inlet.velocity_m_per_s =
        material_inward_velocity[material_inlet_patch];
    material_inlet.thermal_authority =
        hundun::config::InletThermalAuthority::enthalpy;
    material_inlet.enthalpy_J_per_kg = 300000.0;
    material_inlet.density_kg_per_m3 = 1.0;
    material_inlet.scalar_values =
        std::vector<hundun::config::InletScalarValue>{{"alpha", 0.2}};
    auto &material_outlet =
        material_patch_config.boundaries[material_outlet_patch];
    material_outlet.type = hundun::config::BoundaryType::pressure_outlet;
    material_outlet.pressure_perturbation_pa = 0.0;
    auto material_open_boundaries =
        hundun::boundary::BoundaryRegistry::create(material_patch_config,
                                                   open_topology);
    auto material_open_state = make_material_open_state();
    auto material_open_controller =
        hundun::flow::Bdf2RetryController::create(
            material_patch_config.time, material_patch_config.density_model,
            open_topology, open_geometry, mpi, material_open_state);
    const auto material_open_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-r6-material-open-" +
         std::to_string(material_inlet_patch) + "-" +
         std::to_string(mpi.size()));
    if (mpi.rank() == 0)
      std::filesystem::remove_all(material_open_directory);
    mpi.barrier();
    const auto material_open_write = hundun::flow::write_checkpoint_v2(
        mpi, open_decomposition, open_topology, open_geometry,
        material_open_boundaries, material_patch_config, material_open_state,
        material_open_controller.state(), std::nullopt,
        material_open_directory);
    HUNDUN_CHECK(material_open_write.disposition() ==
                 hundun::flow::CheckpointV2Disposition::completed);

    auto changed_material_open_config = material_patch_config;
    const ConfigMutation r5_density_numeric_value =
        [material_inlet_patch](auto &single) {
          single.boundaries[material_inlet_patch].density_kg_per_m3 = 1.1;
        };
    const auto &r6_material_density_all_patches =
        r5_density_numeric_value;
    r6_material_density_all_patches(changed_material_open_config);
    auto changed_material_open_boundaries =
        hundun::boundary::BoundaryRegistry::create(
            changed_material_open_config, open_topology);
    HUNDUN_CHECK(
        independent_resolved_fingerprint(changed_material_open_config) !=
        independent_resolved_fingerprint(material_patch_config));
    HUNDUN_CHECK(
        independent_boundary_fingerprint(changed_material_open_boundaries) !=
        independent_boundary_fingerprint(material_open_boundaries));
    auto changed_material_open_state = make_material_open_destination_state();
    const auto changed_material_open_before =
        CheckpointAccess::snapshot(changed_material_open_state);
    const auto changed_material_open_views =
        CheckpointAccess::density_views(changed_material_open_state);
    const auto changed_material_open_read = hundun::flow::read_checkpoint_v2(
        mpi, open_decomposition, open_topology, open_geometry,
        changed_material_open_boundaries, changed_material_open_config,
        changed_material_open_state, material_open_directory);
    require_constructible_fingerprint_failure(
        mpi, material_open_directory, changed_material_open_before,
        changed_material_open_state, changed_material_open_read);
    for (const auto &view : changed_material_open_views)
      HUNDUN_CHECK(rejects([&] {
        static_cast<void>(view(0, 0, 0, 0));
      }));
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(material_open_directory);
    mpi.barrier();
  }

  auto open_resumed = make_open_state();
  const auto open_read = hundun::flow::read_checkpoint_v2(
      mpi, open_decomposition, open_topology, open_geometry, open_boundaries,
      open_config, open_resumed, open_directory);
  HUNDUN_CHECK(open_read.restored());
  HUNDUN_CHECK(open_read.ideal_gas_closure_state_available());
  HUNDUN_CHECK(open_read.ideal_gas_closure_state().mode ==
               hundun::flow::IdealGasPressureMode::open_fixed);
  HUNDUN_CHECK(!open_read.ideal_gas_closure_state().target_mass_kg.has_value());
  auto open_resumed_controller = hundun::flow::Bdf2RetryController::restore(
      open_config.time, open_config.density_model, open_topology, open_geometry,
      mpi, open_resumed, open_read.time_control_state());
  auto open_resumed_closure = hundun::flow::IdealGasClosure::restore(
      open_topology, open_geometry, open_boundaries, mpi, open_registry,
      open_fields, open_resumed, open_spec,
      open_read.ideal_gas_closure_state());
  auto open_halo = hundun::runtime::HaloExchange::create(
      open_decomposition,
      hundun::runtime::ExchangePlan::create(
          open_decomposition, open_decomposition.local_extent(), 2));
  auto open_resumed_halo = hundun::runtime::HaloExchange::create(
      open_decomposition,
      hundun::runtime::ExchangePlan::create(
          open_decomposition, open_decomposition.local_extent(), 2));
  hundun::execution::CpuReferenceContext open_execution;
  hundun::execution::CpuReferenceContext open_resumed_execution;
  hundun::linear::ConjugateGradientSolver open_momentum(open_execution, mpi),
      open_pressure(open_execution, mpi);
  hundun::linear::ConjugateGradientSolver open_resumed_momentum(
      open_resumed_execution, mpi);
  hundun::linear::ConjugateGradientSolver open_resumed_pressure(
      open_resumed_execution, mpi);
  hundun::linear::JacobiPreconditioner open_mx(open_execution),
      open_my(open_execution), open_mz(open_execution),
      open_pressure_preconditioner(open_execution);
  hundun::linear::JacobiPreconditioner open_resumed_mx(open_resumed_execution),
      open_resumed_my(open_resumed_execution),
      open_resumed_mz(open_resumed_execution),
      open_resumed_pressure_preconditioner(open_resumed_execution);
  hundun::flow::MaterialDensityTransportSpec open_material_spec;
  open_material_spec.enthalpy_density = open_rho_h;
  open_material_spec.scalar_densities = {open_rho_alpha};
  open_material_spec.scalar_diffusivities_kg_per_m_s = {0.0};
  auto open_flow = hundun::flow::FixedStepIdealGasFlow::create(
      open_decomposition, open_topology, open_geometry, open_boundaries, mpi,
      open_execution, open_halo, open_momentum, {&open_mx, &open_my, &open_mz},
      open_pressure, open_pressure_preconditioner, open_registry, open_fields,
      open_material_spec, std::move(open_closure));
  auto open_resumed_flow = hundun::flow::FixedStepIdealGasFlow::create(
      open_decomposition, open_topology, open_geometry, open_boundaries, mpi,
      open_resumed_execution, open_resumed_halo, open_resumed_momentum,
      {&open_resumed_mx, &open_resumed_my, &open_resumed_mz},
      open_resumed_pressure, open_resumed_pressure_preconditioner,
      open_registry, open_fields, open_material_spec,
      std::move(open_resumed_closure));
  const auto open_next =
      open_controller.advance(open_state, open_flow, 0.0, {}, {});
  const auto open_resumed_next = open_resumed_controller.advance(
      open_resumed, open_resumed_flow, 0.0, {}, {});
  HUNDUN_CHECK(open_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(open_resumed_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  for (const auto layer :
       {hundun::flow::FlowLayer::history, hundun::flow::FlowLayer::committed,
        hundun::flow::FlowLayer::trial})
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        open_state.snapshot(layer), open_resumed.snapshot(layer)));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      open_controller.state(), open_resumed_controller.state()));
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      open_flow.closure_state(), open_resumed_flow.closure_state()));

  auto invalid_ideal_values = ideal_values;
  invalid_ideal_values.density.front() *= 1.01;
  auto invalid_ideal_state = hundun::flow::FlowState::create(
      ideal_registry,
      {decomposition.local_extent(), topology.local_face_count()}, ideal_fields,
      metadata);
  invalid_ideal_state.seed_accepted_layers(invalid_ideal_values,
                                           invalid_ideal_values);
  const auto invalid_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-invalid-ideal-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(invalid_directory);
  mpi.barrier();
  const auto invalid_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      invalid_ideal_state, ideal_controller_state, persisted_closure,
      invalid_directory);
  HUNDUN_CHECK(invalid_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::failed);
  HUNDUN_CHECK(invalid_write.reason() ==
               hundun::flow::CheckpointV2FailureReason::state);
  HUNDUN_CHECK(!std::filesystem::exists(invalid_directory));

  const auto common_mismatch_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-common-mismatch-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(common_mismatch_directory);
  mpi.barrier();
  if (mpi.size() > 1) {
    auto mismatched_config = config;
    if (mpi.rank() == mpi.size() - 1)
      mismatched_config.physics.dynamic_viscosity_pa_s *= 2.0;
    const auto source_before = CheckpointAccess::snapshot(source);
    const auto mismatch = hundun::flow::write_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, mismatched_config,
        source, controller_state, std::nullopt, common_mismatch_directory);
    HUNDUN_CHECK(mismatch.disposition() ==
                 hundun::flow::CheckpointV2Disposition::failed);
    HUNDUN_CHECK(mismatch.reason() ==
                 hundun::flow::CheckpointV2FailureReason::invalid_input);
    HUNDUN_CHECK(mismatch.phase() ==
                 hundun::flow::CheckpointV2Phase::preflight);
    HUNDUN_CHECK(mismatch.lowest_failing_rank() == mpi.size() - 1);
    HUNDUN_CHECK(mismatch.step() == metadata.step);
    HUNDUN_CHECK(mismatch.time_s() == metadata.time_s);
    HUNDUN_CHECK(!std::filesystem::exists(common_mismatch_directory));
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_deep_snapshot_equal(
        source_before, CheckpointAccess::snapshot(source)));
  }

  const auto nonroot_invalid_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-nonroot-invalid-time-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(nonroot_invalid_directory);
  mpi.barrier();
  if (mpi.size() > 1) {
    auto invalid_time = controller_state;
    if (mpi.rank() == 1)
      invalid_time.proposed_next_dt_s = -1.0;
    const auto nonroot_invalid = hundun::flow::write_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config, source,
        invalid_time, std::nullopt, nonroot_invalid_directory);
    HUNDUN_CHECK(nonroot_invalid.disposition() ==
                 hundun::flow::CheckpointV2Disposition::failed);
    HUNDUN_CHECK(nonroot_invalid.reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(nonroot_invalid.phase() ==
                 hundun::flow::CheckpointV2Phase::preflight);
    HUNDUN_CHECK(nonroot_invalid.lowest_failing_rank() == 1);
    HUNDUN_CHECK(nonroot_invalid.step() == metadata.step);
    HUNDUN_CHECK(nonroot_invalid.time_s() == metadata.time_s);
    HUNDUN_CHECK(!std::filesystem::exists(nonroot_invalid_directory));
  }

  using ReadFault = hundun::runtime::checkpoint_v2::test::ExactReadFault;
  struct WritePhaseCase final {
    std::uint32_t calls_before;
    int failing_rank;
    hundun::flow::CheckpointV2Phase phase;
  };
  const std::array write_phase_cases{
      WritePhaseCase{0U, mpi.size() - 1,
                     hundun::flow::CheckpointV2Phase::rank_temporary_file},
      WritePhaseCase{1U, 0, hundun::flow::CheckpointV2Phase::manifest},
      WritePhaseCase{2U, 0, hundun::flow::CheckpointV2Phase::completed_marker}};
  for (std::size_t index = 0; index < write_phase_cases.size(); ++index) {
    const auto failure_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-write-phase-" + std::to_string(index) + "-" +
         std::to_string(mpi.size()));
    if (mpi.rank() == 0)
      std::filesystem::remove_all(failure_directory);
    mpi.barrier();
    if (mpi.rank() == write_phase_cases[index].failing_rank)
      hundun::runtime::checkpoint_v2::test::set_exact_read_fault(
          ReadFault::truncate_after_size,
          write_phase_cases[index].calls_before);
    const auto failed_phase = hundun::flow::write_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config, source,
        controller_state, std::nullopt, failure_directory);
    HUNDUN_CHECK(failed_phase.disposition() ==
                 hundun::flow::CheckpointV2Disposition::failed);
    HUNDUN_CHECK(failed_phase.reason() ==
                 hundun::flow::CheckpointV2FailureReason::file_integrity);
    HUNDUN_CHECK(failed_phase.phase() == write_phase_cases[index].phase);
    HUNDUN_CHECK(failed_phase.lowest_failing_rank() ==
                 write_phase_cases[index].failing_rank);
    HUNDUN_CHECK(failed_phase.exact_size_and_eof_status() ==
                 hundun::flow::CheckpointV2CheckStatus::failed);
    HUNDUN_CHECK(!std::filesystem::exists(failure_directory / "COMPLETED"));
    require_consistent_product_report(mpi, failed_phase);
    auto expected_integrity =
        expected_write_failure(
            write_phase_cases[index].phase,
            write_phase_cases[index].failing_rank,
            index == 0U
                ? 0U
                : index == 1U
                      ? static_cast<std::uint64_t>(mpi.size())
                      : static_cast<std::uint64_t>(mpi.size()) + 1U,
            index == 0U
                ? static_cast<std::uint64_t>(mpi.size() - 1)
                : index == 1U
                      ? static_cast<std::uint64_t>(mpi.size())
                      : static_cast<std::uint64_t>(mpi.size()) + 1U,
            index == 0U ? 22U : index == 1U ? 30U : 33U,
            index == 0U &&
                    mpi.rank() == write_phase_cases[index].failing_rank
                ? hundun::flow::CheckpointV2CheckStatus::failed
                : hundun::flow::CheckpointV2CheckStatus::passed,
            hundun::flow::CheckpointV2CheckStatus::failed,
            index == 0U
                ? hundun::flow::CheckpointV2CheckStatus::not_checked
                : hundun::flow::CheckpointV2CheckStatus::failed);
    expected_integrity.values.reason =
        hundun::flow::CheckpointV2FailureReason::file_integrity;
    if (index >= 1U)
      expected_integrity.values.manifest_crc =
          index == 1U ? hundun::flow::CheckpointV2CheckStatus::failed
                      : hundun::flow::CheckpointV2CheckStatus::passed;
    if (index == 2U) {
      expected_integrity.values.global_logical_bytes =
          expected_global_logical;
      expected_integrity.values.global_actual_bytes =
          expected_rank_actual_sum + expected_manifest_size + 40U;
      expected_integrity.values.manifest_crc64 =
          independent_crc64(manifest_bytes_for_report);
      expected_integrity.values.fingerprint =
          hundun::flow::CheckpointV2CheckStatus::failed;
    }
    expected_integrity.semantic_fingerprint =
        independent_report_fingerprint(expected_integrity.values);
    require_exact_product_report(failed_phase, expected_integrity);
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(failure_directory);
    mpi.barrier();
  }

  auto constant_config = config;
  constant_config.density_model = hundun::config::DensityModel::constant;
  auto constant_boundaries =
      hundun::boundary::BoundaryRegistry::create(constant_config, topology);
  auto constant_values = committed;
  constant_values.density.assign(cells, 1.0);
  constant_values.velocity.assign(cells * 3U, 0.0);
  constant_values.mechanical_pressure.assign(cells, 0.0);
  constant_values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  constant_values.face_mass_flux.assign(topology.local_face_count(), 0.0);
  constant_values.transported_cell_fields = {
      std::vector<double>(cells, 300000.0), std::vector<double>(cells, 0.2)};
  auto uninterrupted = make_state();
  uninterrupted.seed_accepted_layers(constant_values, constant_values);
  auto uninterrupted_controller = hundun::flow::Bdf2RetryController::create(
      constant_config.time, constant_config.density_model, topology, geometry,
      mpi, uninterrupted);
  const auto continuation_controller_state = uninterrupted_controller.state();
  const auto continuation_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-continuation-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(continuation_directory);
  mpi.barrier();
  const auto continuation_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, constant_boundaries,
      constant_config, uninterrupted, continuation_controller_state,
      std::nullopt, continuation_directory);
  HUNDUN_CHECK(continuation_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);

  const auto bytes_u8 = [](std::uint8_t value) { return Bytes{value}; };
  const auto bytes_u32 = [](std::uint32_t value) {
    Bytes bytes;
    append_u32(bytes, value);
    return bytes;
  };
  const auto bytes_u64 = [](std::uint64_t value) {
    Bytes bytes;
    append_u64(bytes, value);
    return bytes;
  };
  const auto bytes_f64 = [](double value) {
    Bytes bytes;
    append_f64(bytes, value);
    return bytes;
  };
  struct TimeStateMutation final {
    const char *name;
    std::size_t offset;
    Bytes replacement;
    hundun::flow::CheckpointV2FailureReason reason;
    hundun::flow::CheckpointV2Phase phase;
  };
  const std::array time_state_mutations{
      TimeStateMutation{
          "schema", 37U,
          bytes_u32(continuation_controller_state.schema_version + 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{
          "accepted-step", 41U,
          bytes_u64(continuation_controller_state.accepted_step + 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"proposed-dt", 49U, bytes_f64(-1.0),
                        hundun::flow::CheckpointV2FailureReason::state,
                        hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"last-dt", 57U, bytes_f64(-1.0),
                        hundun::flow::CheckpointV2FailureReason::state,
                        hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"order", 65U, bytes_u8(2U),
                        hundun::flow::CheckpointV2FailureReason::file_integrity,
                        hundun::flow::CheckpointV2Phase::manifest_read},
      TimeStateMutation{
          "history", 66U,
          bytes_u8(continuation_controller_state.history_ready ? 0U : 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{
          "linear-half", 67U,
          bytes_u8(continuation_controller_state
                           .last_all_linear_solves_within_half_limit
                       ? 0U
                       : 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"convective-rate", 68U,
                        bytes_u64(UINT64_C(0x7ff8000000000001)),
                        hundun::flow::CheckpointV2FailureReason::state,
                        hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"diffusive-rate", 76U,
                        bytes_u64(UINT64_C(0x7ff8000000000001)),
                        hundun::flow::CheckpointV2FailureReason::state,
                        hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{
          "metrics", 84U,
          bytes_u8(
              continuation_controller_state.last_stability_metrics_available
                  ? 0U
                  : 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{
          "retry-count", 85U,
          bytes_u32(static_cast<std::uint32_t>(config.time.max_retries + 1)),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{"revision", 89U,
                        bytes_u64(continuation_controller_state.revision + 1U),
                        hundun::flow::CheckpointV2FailureReason::state,
                        hundun::flow::CheckpointV2Phase::restore_prepare},
      TimeStateMutation{
          "seal", 97U, bytes_u64(continuation_controller_state.state_seal ^ 1U),
          hundun::flow::CheckpointV2FailureReason::state,
          hundun::flow::CheckpointV2Phase::restore_prepare}};
  for (const auto &mutation : time_state_mutations) {
    const auto invalid_time_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-time-state-" + std::string(mutation.name) + "-" +
         std::to_string(mpi.size()));
    if (mpi.rank() == 0) {
      std::filesystem::remove_all(invalid_time_directory);
      std::filesystem::copy(continuation_directory, invalid_time_directory,
                            std::filesystem::copy_options::recursive);
    }
    mpi.barrier();
    rebuild_checkpoint_global_bytes(mpi, invalid_time_directory,
                                    mutation.offset, mutation.replacement);
    auto invalid_time_destination = make_destination_state();
    invalid_time_destination.seed_accepted_layers(different, different);
    const auto before = CheckpointAccess::snapshot(invalid_time_destination);
    const auto rejected = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, constant_boundaries,
        constant_config, invalid_time_destination, invalid_time_directory);
    HUNDUN_CHECK(!rejected.restored());
    HUNDUN_CHECK(rejected.report().reason() == mutation.reason);
    HUNDUN_CHECK(rejected.report().phase() == mutation.phase);
    require_destination_report_authority(mpi, before, rejected.report());
    if (std::string_view(mutation.name) == "schema")
      require_exact_product_report(
          rejected.report(),
          expected_late_read_report(invalid_time_directory, before,
                                    LateReadEvidence::global_state, 0));
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        before, CheckpointAccess::snapshot(invalid_time_destination)));
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(invalid_time_directory);
  }

  struct ReadPhaseCase final {
    std::uint32_t calls_before;
    hundun::flow::CheckpointV2Phase phase;
  };
  constexpr std::array read_phase_cases{
      ReadPhaseCase{0U, hundun::flow::CheckpointV2Phase::marker_read},
      ReadPhaseCase{1U, hundun::flow::CheckpointV2Phase::manifest_read},
      ReadPhaseCase{2U, hundun::flow::CheckpointV2Phase::rank_read}};
  for (const auto &phase_case : read_phase_cases) {
    auto failed_destination = make_destination_state();
    failed_destination.seed_accepted_layers(different, different);
    const auto before_failure = CheckpointAccess::snapshot(failed_destination);
    const int failing_rank = mpi.size() - 1;
    if (mpi.rank() == failing_rank)
      hundun::runtime::checkpoint_v2::test::set_exact_read_fault(
          ReadFault::read_failure, phase_case.calls_before);
    const auto failed_phase = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, constant_boundaries,
        constant_config, failed_destination, continuation_directory);
    HUNDUN_CHECK(!failed_phase.restored());
    HUNDUN_CHECK(failed_phase.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::filesystem);
    HUNDUN_CHECK(failed_phase.report().phase() == phase_case.phase);
    HUNDUN_CHECK(failed_phase.report().lowest_failing_rank() == failing_rank);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        before_failure, CheckpointAccess::snapshot(failed_destination)));
    require_destination_report_authority(mpi, before_failure,
                                         failed_phase.report());
    auto expected_read_fault =
        expected_read_file_failure(phase_case.phase, failing_rank);
    if (phase_case.phase == hundun::flow::CheckpointV2Phase::rank_read) {
      const auto continuation_manifest =
          read_file_bytes(continuation_directory / "manifest.v2.bin");
      expected_read_fault.values.manifest_crc64 =
          independent_crc64(continuation_manifest);
      if (mpi.rank() != failing_rank) {
        const auto continuation_rank =
            read_file_bytes(rank_file_path(continuation_directory, mpi.rank()));
        expected_read_fault.values.local_actual_bytes =
            continuation_rank.size();
        expected_read_fault.values.local_crc64 =
            independent_crc64(continuation_rank);
      }
      expected_read_fault.semantic_fingerprint =
          independent_report_fingerprint(expected_read_fault.values);
    }
    require_exact_product_report(failed_phase.report(),
                                 expected_read_fault);
  }

  constexpr std::array<std::uint64_t, 4> invalid_density_bits{
      UINT64_C(0xbff0000000000000), UINT64_C(0x0000000000000000),
      UINT64_C(0x7ff8000000000001), UINT64_C(0x4000000000000000)};
  for (std::size_t invalid_index = 0;
       invalid_index < invalid_density_bits.size(); ++invalid_index) {
    const auto invalid_state_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-valid-crc-invalid-state-" +
         std::to_string(invalid_index) + "-" + std::to_string(mpi.size()));
    if (mpi.rank() == 0) {
      std::filesystem::remove_all(invalid_state_directory);
      std::filesystem::copy(continuation_directory, invalid_state_directory,
                            std::filesystem::copy_options::recursive);
    }
    mpi.barrier();
    rebuild_checkpoint_rank_value(mpi, invalid_state_directory,
                                  topology.owned_cell_count(),
                                  topology.local_face_count(), 2U, 0U,
                                  invalid_density_bits[invalid_index]);
    auto invalid_destination = make_state();
    invalid_destination.seed_accepted_layers(different, different);
    const auto before_invalid = CheckpointAccess::snapshot(invalid_destination);
    const auto invalid_read = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, constant_boundaries,
        constant_config, invalid_destination, invalid_state_directory);
    HUNDUN_CHECK(!invalid_read.restored());
    HUNDUN_CHECK(invalid_read.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(invalid_read.report().phase() ==
                 hundun::flow::CheckpointV2Phase::restore_prepare);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        before_invalid, CheckpointAccess::snapshot(invalid_destination)));
    require_consistent_product_report(mpi, invalid_read.report());
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(invalid_state_directory);
    mpi.barrier();
  }

  const auto invalid_generic_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-valid-crc-invalid-generic-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(invalid_generic_directory);
    std::filesystem::copy(continuation_directory, invalid_generic_directory,
                          std::filesystem::copy_options::recursive);
  }
  mpi.barrier();
  rebuild_checkpoint_rank_value(
      mpi, invalid_generic_directory, topology.owned_cell_count(),
      topology.local_face_count(), 2U, 6U, UINT64_C(0x7ff8000000000001));
  auto invalid_generic_destination = make_destination_state();
  invalid_generic_destination.seed_accepted_layers(different, different);
  const auto before_invalid_generic =
      CheckpointAccess::snapshot(invalid_generic_destination);
  const auto invalid_generic_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, constant_boundaries,
      constant_config, invalid_generic_destination, invalid_generic_directory);
  HUNDUN_CHECK(!invalid_generic_read.restored());
  HUNDUN_CHECK(invalid_generic_read.report().reason() ==
               hundun::flow::CheckpointV2FailureReason::state);
  HUNDUN_CHECK(invalid_generic_read.report().phase() ==
               hundun::flow::CheckpointV2Phase::restore_prepare);
  HUNDUN_CHECK(invalid_generic_read.report().lowest_failing_rank() ==
               mpi.size() - 1);
  HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
      before_invalid_generic,
      CheckpointAccess::snapshot(invalid_generic_destination)));
  require_destination_report_authority(mpi, before_invalid_generic,
                                       invalid_generic_read.report());
  require_exact_product_report(
      invalid_generic_read.report(),
      expected_late_read_report(invalid_generic_directory,
                                before_invalid_generic,
                                LateReadEvidence::physical_state,
                                mpi.size() - 1));
  mpi.barrier();
  if (mpi.rank() == 0)
    std::filesystem::remove_all(invalid_generic_directory);
  mpi.barrier();

  struct IdealStateMutation final {
    const char *name{};
    std::size_t record{};
    std::uint64_t bits{};
  };
  const std::array ideal_state_mutations{
      IdealStateMutation{"enthalpy", 5U, UINT64_C(0x0000000000000000)},
      IdealStateMutation{
          "eos", 5U,
          bits_of(ideal_values.transported_cell_fields[0][0] * 1.01)},
      IdealStateMutation{"closed-mass", 0U,
                         bits_of(ideal_values.density[0] * 1.01)}};
  for (const auto &mutation : ideal_state_mutations) {
    const auto invalid_ideal_directory =
        std::filesystem::temp_directory_path() /
        ("hundun-task23-valid-crc-invalid-" + std::string(mutation.name) + "-" +
         std::to_string(mpi.size()));
    if (mpi.rank() == 0) {
      std::filesystem::remove_all(invalid_ideal_directory);
      std::filesystem::copy(ideal_directory, invalid_ideal_directory,
                            std::filesystem::copy_options::recursive);
    }
    mpi.barrier();
    rebuild_checkpoint_rank_value(
        mpi, invalid_ideal_directory, topology.owned_cell_count(),
        topology.local_face_count(), 1U, mutation.record, mutation.bits);
    auto invalid_ideal_destination = hundun::flow::FlowState::create(
        ideal_registry,
        {decomposition.local_extent(), topology.local_face_count()},
        ideal_fields, metadata);
    invalid_ideal_destination.seed_accepted_layers(ideal_values, ideal_values);
    const auto before_invalid_ideal =
        CheckpointAccess::snapshot(invalid_ideal_destination);
    const auto invalid_ideal_read = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
        invalid_ideal_destination, invalid_ideal_directory);
    HUNDUN_CHECK(!invalid_ideal_read.restored());
    HUNDUN_CHECK(invalid_ideal_read.report().reason() ==
                 hundun::flow::CheckpointV2FailureReason::state);
    HUNDUN_CHECK(invalid_ideal_read.report().phase() ==
                 hundun::flow::CheckpointV2Phase::restore_prepare);
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
        before_invalid_ideal,
        CheckpointAccess::snapshot(invalid_ideal_destination)));
    require_consistent_product_report(mpi, invalid_ideal_read.report());
    mpi.barrier();
    if (mpi.rank() == 0)
      std::filesystem::remove_all(invalid_ideal_directory);
    mpi.barrier();
  }

  auto resumed = make_state();
  resumed.seed_accepted_layers(different, different);
  const auto continuation_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, constant_boundaries,
      constant_config, resumed, continuation_directory);
  HUNDUN_CHECK(continuation_read.restored());
  auto resumed_controller = hundun::flow::Bdf2RetryController::restore(
      constant_config.time, constant_config.density_model, topology, geometry,
      mpi, resumed, continuation_read.time_control_state());

  const auto local_extent = decomposition.local_extent();
  auto uninterrupted_halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local_extent, 2));
  auto resumed_halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local_extent, 2));
  hundun::execution::CpuReferenceContext uninterrupted_execution;
  hundun::execution::CpuReferenceContext resumed_execution;
  hundun::linear::ConjugateGradientSolver uninterrupted_momentum(
      uninterrupted_execution, mpi);
  hundun::linear::ConjugateGradientSolver uninterrupted_pressure(
      uninterrupted_execution, mpi);
  hundun::linear::ConjugateGradientSolver resumed_momentum(resumed_execution,
                                                           mpi);
  hundun::linear::ConjugateGradientSolver resumed_pressure(resumed_execution,
                                                           mpi);
  hundun::linear::JacobiPreconditioner uninterrupted_mx(
      uninterrupted_execution),
      uninterrupted_my(uninterrupted_execution),
      uninterrupted_mz(uninterrupted_execution),
      uninterrupted_pressure_preconditioner(uninterrupted_execution);
  hundun::linear::JacobiPreconditioner resumed_mx(resumed_execution),
      resumed_my(resumed_execution), resumed_mz(resumed_execution),
      resumed_pressure_preconditioner(resumed_execution);
  const std::vector<hundun::flow::ConstantDensityTransportSpec> transported{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::enthalpy(), 0.0},
      {fields.transported_cell_fields[1],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}};
  auto uninterrupted_flow = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, constant_boundaries, mpi,
      uninterrupted_execution, uninterrupted_halo, uninterrupted_momentum,
      {&uninterrupted_mx, &uninterrupted_my, &uninterrupted_mz},
      uninterrupted_pressure, uninterrupted_pressure_preconditioner,
      transported);
  auto resumed_flow = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, constant_boundaries, mpi,
      resumed_execution, resumed_halo, resumed_momentum,
      {&resumed_mx, &resumed_my, &resumed_mz}, resumed_pressure,
      resumed_pressure_preconditioner, transported);
  const auto uninterrupted_step = uninterrupted_controller.advance(
      uninterrupted, uninterrupted_flow, 1.0, 0.0, {}, {});
  const auto resumed_step =
      resumed_controller.advance(resumed, resumed_flow, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(uninterrupted_step.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(resumed_step.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  for (const auto layer :
       {hundun::flow::FlowLayer::history, hundun::flow::FlowLayer::committed,
        hundun::flow::FlowLayer::trial}) {
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        uninterrupted.snapshot(layer), resumed.snapshot(layer)));
  }
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      uninterrupted.metadata(), resumed.metadata()));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      uninterrupted_controller.state(), resumed_controller.state()));

  auto material_uninterrupted = make_state();
  material_uninterrupted.seed_accepted_layers(constant_values, constant_values);
  auto material_controller = hundun::flow::Bdf2RetryController::create(
      config.time, config.density_model, topology, geometry, mpi,
      material_uninterrupted);
  const auto material_continuation_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-material-continuation-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(material_continuation_directory);
  mpi.barrier();
  const auto material_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config,
      material_uninterrupted, material_controller.state(), std::nullopt,
      material_continuation_directory);
  HUNDUN_CHECK(material_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  auto material_resumed = make_state();
  material_resumed.seed_accepted_layers(different, different);
  const auto material_read = hundun::flow::read_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config,
      material_resumed, material_continuation_directory);
  HUNDUN_CHECK(material_read.restored());
  auto material_resumed_controller = hundun::flow::Bdf2RetryController::restore(
      config.time, config.density_model, topology, geometry, mpi,
      material_resumed, material_read.time_control_state());
  hundun::flow::MaterialDensityTransportSpec material_spec;
  material_spec.enthalpy_density = fields.transported_cell_fields.front();
  material_spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  material_spec.scalar_densities = {fields.transported_cell_fields[1]};
  material_spec.scalar_diffusivities_kg_per_m_s = {0.0};
  auto material_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi,
      uninterrupted_execution, uninterrupted_halo, uninterrupted_momentum,
      {&uninterrupted_mx, &uninterrupted_my, &uninterrupted_mz},
      uninterrupted_pressure, uninterrupted_pressure_preconditioner, registry,
      fields, material_spec);
  auto material_resumed_flow =
      hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, resumed_execution,
          resumed_halo, resumed_momentum,
          {&resumed_mx, &resumed_my, &resumed_mz}, resumed_pressure,
          resumed_pressure_preconditioner, registry, fields, material_spec);
  const auto material_next = material_controller.advance(
      material_uninterrupted, material_flow, 0.0, {}, {});
  const auto material_resumed_next = material_resumed_controller.advance(
      material_resumed, material_resumed_flow, 0.0, {}, {});
  HUNDUN_CHECK(material_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(material_resumed_next.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  for (const auto layer :
       {hundun::flow::FlowLayer::history, hundun::flow::FlowLayer::committed,
        hundun::flow::FlowLayer::trial}) {
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        material_uninterrupted.snapshot(layer),
        material_resumed.snapshot(layer)));
  }
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      material_uninterrupted.metadata(), material_resumed.metadata()));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      material_controller.state(), material_resumed_controller.state()));

  mpi.barrier();
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(directory);
    std::filesystem::remove_all(ideal_directory);
    std::filesystem::remove_all(ideal_max_revision_directory);
    std::filesystem::remove_all(open_directory);
    std::filesystem::remove_all(invalid_directory);
    std::filesystem::remove_all(common_mismatch_directory);
    std::filesystem::remove_all(nonroot_invalid_directory);
    std::filesystem::remove_all(continuation_directory);
    std::filesystem::remove_all(material_continuation_directory);
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::string mode(argv[1]);
  if (mode != "fast" && mode != "acceptance")
    return 2;
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi, mode == "acceptance");
  });
}

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

std::uint8_t code(config::SimulationType item) {
  switch (item) {
  case config::SimulationType::passive_scalar: return 0U;
  case config::SimulationType::variable_density_flow: return 1U;
  }
  throw runtime::Error("Checkpoint v2 simulation type is invalid");
}
std::uint8_t code(config::DensityModel item) {
  switch (item) {
  case config::DensityModel::constant: return 0U;
  case config::DensityModel::material: return 1U;
  case config::DensityModel::ideal_gas: return 2U;
  }
  throw runtime::Error("Checkpoint v2 density model is invalid");
}
std::uint8_t code(config::MeshMapping item) {
  switch (item) {
  case config::MeshMapping::uniform_box: return 0U;
  case config::MeshMapping::analytic_warped_box: return 1U;
  }
  throw runtime::Error("Checkpoint v2 mesh mapping is invalid");
}
std::uint8_t code(config::TimeMode item) {
  switch (item) {
  case config::TimeMode::fixed: return 0U;
  case config::TimeMode::adaptive: return 1U;
  }
  throw runtime::Error("Checkpoint v2 time mode is invalid");
}
std::uint8_t code(config::PatchName item) {
  switch (item) {
  case config::PatchName::x_min: return 0U;
  case config::PatchName::x_max: return 1U;
  case config::PatchName::y_min: return 2U;
  case config::PatchName::y_max: return 3U;
  case config::PatchName::z_min: return 4U;
  case config::PatchName::z_max: return 5U;
  }
  throw runtime::Error("Checkpoint v2 patch code is invalid");
}
std::uint8_t code(config::BoundaryType item) {
  switch (item) {
  case config::BoundaryType::periodic: return 0U;
  case config::BoundaryType::no_slip_wall: return 1U;
  case config::BoundaryType::symmetry: return 2U;
  case config::BoundaryType::velocity_inlet: return 3U;
  case config::BoundaryType::pressure_outlet: return 4U;
  }
  throw runtime::Error("Checkpoint v2 boundary type is invalid");
}
std::uint8_t code(config::InletThermalAuthority item) {
  switch (item) {
  case config::InletThermalAuthority::temperature: return 0U;
  case config::InletThermalAuthority::enthalpy: return 1U;
  }
  throw runtime::Error("Checkpoint v2 inlet authority is invalid");
}
std::uint8_t code(runtime::FunctionSpace item) {
  switch (item) {
  case runtime::FunctionSpace::cell_average: return 0U;
  case runtime::FunctionSpace::face_value: return 1U;
  case runtime::FunctionSpace::vertex_value: return 2U;
  case runtime::FunctionSpace::element_dof: return 3U;
  case runtime::FunctionSpace::quadrature_point: return 4U;
  case runtime::FunctionSpace::particle: return 5U;
  }
  throw runtime::Error("Checkpoint v2 function space is invalid");
}
std::uint8_t code(runtime::ScalarType item) {
  switch (item) {
  case runtime::ScalarType::float64: return 0U;
  case runtime::ScalarType::int32: return 1U;
  case runtime::ScalarType::uint8: return 2U;
  }
  throw runtime::Error("Checkpoint v2 scalar type is invalid");
}
std::uint8_t code(runtime::RestartPolicy item) {
  switch (item) {
  case runtime::RestartPolicy::persistent: return 0U;
  case runtime::RestartPolicy::transient: return 1U;
  }
  throw runtime::Error("Checkpoint v2 restart policy is invalid");
}
std::uint8_t code(runtime::OutputPolicy item) {
  switch (item) {
  case runtime::OutputPolicy::never: return 0U;
  case runtime::OutputPolicy::selected: return 1U;
  case runtime::OutputPolicy::always: return 2U;
  }
  throw runtime::Error("Checkpoint v2 output policy is invalid");
}
std::uint8_t code(boundary::BoundaryKind item) {
  switch (item) {
  case boundary::BoundaryKind::periodic: return 0U;
  case boundary::BoundaryKind::no_slip_wall: return 1U;
  case boundary::BoundaryKind::symmetry: return 2U;
  case boundary::BoundaryKind::velocity_inlet: return 3U;
  case boundary::BoundaryKind::pressure_outlet: return 4U;
  }
  throw runtime::Error("Checkpoint v2 resolved boundary kind is invalid");
}
std::uint8_t code(boundary::VelocityRule item) {
  switch (item) {
  case boundary::VelocityRule::periodic_pair: return 0U;
  case boundary::VelocityRule::prescribed_zero: return 1U;
  case boundary::VelocityRule::reflect_normal_copy_tangential: return 2U;
  case boundary::VelocityRule::prescribed_inlet: return 3U;
  case boundary::VelocityRule::pure_outflow: return 4U;
  }
  throw runtime::Error("Checkpoint v2 velocity rule is invalid");
}
std::uint8_t code(boundary::PressureRule item) {
  switch (item) {
  case boundary::PressureRule::periodic_pair: return 0U;
  case boundary::PressureRule::zero_normal_gradient: return 1U;
  case boundary::PressureRule::prescribed_value: return 2U;
  }
  throw runtime::Error("Checkpoint v2 pressure rule is invalid");
}
std::uint8_t code(boundary::TransportRule item) {
  switch (item) {
  case boundary::TransportRule::periodic_pair: return 0U;
  case boundary::TransportRule::copy_interior: return 1U;
  case boundary::TransportRule::zero_normal_diffusive_flux: return 2U;
  case boundary::TransportRule::prescribed_value: return 3U;
  case boundary::TransportRule::pure_outflow: return 4U;
  }
  throw runtime::Error("Checkpoint v2 transport rule is invalid");
}
std::uint8_t code(boundary::MassFluxRule item) {
  switch (item) {
  case boundary::MassFluxRule::periodic_pair: return 0U;
  case boundary::MassFluxRule::identically_zero: return 1U;
  case boundary::MassFluxRule::prescribed_inlet_state: return 2U;
  case boundary::MassFluxRule::outflow_only: return 3U;
  }
  throw runtime::Error("Checkpoint v2 mass-flux rule is invalid");
}
std::uint8_t code(MomentumTimeOrder item) {
  switch (item) {
  case MomentumTimeOrder::backward_euler: return 0U;
  case MomentumTimeOrder::bdf2: return 1U;
  }
  throw runtime::Error("Checkpoint v2 momentum order is invalid");
}
std::uint8_t code(IdealGasPressureMode item) {
  switch (item) {
  case IdealGasPressureMode::closed_dynamic: return 0U;
  case IdealGasPressureMode::open_fixed: return 1U;
  }
  throw runtime::Error("Checkpoint v2 pressure mode is invalid");
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
void append(Encoder &encoder, runtime::Int3 item) {
  encoder.i32(item.x);
  encoder.i32(item.y);
  encoder.i32(item.z);
}
void append(Encoder &encoder, runtime::Box3 item) {
  append(encoder, item.begin);
  append(encoder, item.end);
}
void append(Encoder &encoder, runtime::Real3 item) {
  encoder.f64(item.x);
  encoder.f64(item.y);
  encoder.f64(item.z);
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
    encoder.u8(code(config.simulation_type));
    encoder.u8(code(config.density_model));
    encoder.boolean(config.resources.expected_ranks.has_value());
    if (config.resources.expected_ranks)
      encoder.i32(*config.resources.expected_ranks);
    encoder.boolean(config.resources.process_grid.has_value());
    if (config.resources.process_grid) {
      encoder.i32(config.resources.process_grid->x);
      encoder.i32(config.resources.process_grid->y);
      encoder.i32(config.resources.process_grid->z);
    }
    append(encoder, config.mesh.cells);
    append(encoder, config.mesh.origin_m);
    append(encoder, config.mesh.length_m);
    encoder.u8(code(config.mesh.mapping));
    append_optional(encoder, config.mesh.warp_amplitude);
    encoder.u8(code(config.time.mode));
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
    for (const auto &boundary : config.boundaries) {
      encoder.u8(code(boundary.patch));
      encoder.u8(code(boundary.type));
      append_optional(encoder, boundary.velocity_m_per_s);
      encoder.boolean(boundary.thermal_authority.has_value());
      if (boundary.thermal_authority)
        encoder.u8(code(*boundary.thermal_authority));
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

std::uint64_t local_topology_fingerprint(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, int rank, int rank_count) {
  return make_fingerprint("hundun.checkpoint-v2.topology-local.v1",
                          [&](Encoder &encoder) {
    encoder.i32(rank);
    encoder.i32(rank_count);
    append(encoder, decomposition.process_grid());
    append(encoder, topology.global_extent());
    append(encoder, topology.owned_global_box());
    encoder.u64(topology.global_cell_count());
    encoder.u64(topology.global_face_count());
    for (const auto axis : {mesh::FaceAxis::x, mesh::FaceAxis::y,
                            mesh::FaceAxis::z})
      encoder.u64(topology.global_face_count(axis));
    encoder.u64(topology.owned_cell_count());
    encoder.u64(topology.ghost_cell_count());
    encoder.u64(topology.local_cell_count());
    encoder.u64(topology.owned_face_count());
    encoder.u64(topology.ghost_face_count());
    encoder.u64(topology.local_face_count());
    for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
      encoder.u64(topology.global_cell_id(cell));
      append(encoder, topology.global_cell(cell));
      encoder.u8(topology.cell_ownership(cell) == mesh::EntityOwnership::owned
                     ? 0U : 1U);
    }
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      encoder.u64(topology.global_face_id(face));
      const auto logical = topology.logical_face(face);
      encoder.u8(logical.axis == mesh::FaceAxis::x
                     ? 0U : logical.axis == mesh::FaceAxis::y ? 1U : 2U);
      append(encoder, logical.coordinate);
      encoder.u8(topology.face_ownership(face) == mesh::EntityOwnership::owned
                     ? 0U : 1U);
      encoder.u64(topology.global_cell_id(topology.owner(face)));
      const auto neighbour = topology.neighbour(face);
      encoder.boolean(neighbour.has_value());
      if (neighbour)
        encoder.u64(topology.global_cell_id(*neighbour));
      const auto patch = topology.patch_id(face);
      encoder.boolean(patch.has_value());
      if (patch)
        encoder.u32(*patch);
      const auto pair = topology.periodic_pair(face);
      encoder.boolean(pair.has_value());
      if (pair)
        encoder.u64(*pair);
    }
    for (const auto &patch : topology.patches()) {
      encoder.u32(patch.stable_id());
      encoder.string(std::string(patch.name()));
      encoder.u8(patch.pairing_kind() == mesh::PatchPairingKind::none ? 0U
                                                                      : 1U);
      encoder.boolean(patch.paired_patch_id().has_value());
      if (patch.paired_patch_id())
        encoder.u32(*patch.paired_patch_id());
      encoder.u64(patch.local_faces().size());
      for (const auto face : patch.local_faces())
        encoder.u64(face);
    }
  });
}

std::uint64_t local_geometry_fingerprint(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const config::FlowCaseConfig &config, int rank) {
  return make_fingerprint("hundun.checkpoint-v2.geometry-local.v1",
                          [&](Encoder &encoder) {
    encoder.i32(rank);
    encoder.u8(geometry.mapping_kind() == mesh::MappingKind::uniform_box
                   ? 0U : 1U);
    append(encoder, geometry.global_extent());
    append(encoder, geometry.owned_global_box());
    append(encoder, geometry.origin_m());
    append(encoder, geometry.length_m());
    append_optional(encoder, geometry.uniform_spacing_m());
    append_optional(encoder, config.mesh.warp_amplitude);
    const auto box = topology.owned_global_box();
    if (box.end.x < box.begin.x || box.end.y < box.begin.y ||
        box.end.z < box.begin.z)
      throw runtime::Error("Checkpoint v2 geometry owned box is invalid");
    auto vertex_count = runtime::checkpoint_v2::checked_product(
        static_cast<std::size_t>(box.end.x - box.begin.x + 1),
        static_cast<std::size_t>(box.end.y - box.begin.y + 1));
    vertex_count = runtime::checkpoint_v2::checked_product(
        vertex_count,
        static_cast<std::size_t>(box.end.z - box.begin.z + 1));
    encoder.u64(vertex_count);
    for (int k = box.begin.z; k <= box.end.z; ++k)
      for (int j = box.begin.y; j <= box.end.y; ++j)
        for (int i = box.begin.x; i <= box.end.x; ++i) {
          const runtime::Int3 coordinate{i, j, k};
          append(encoder, coordinate);
          append(encoder, geometry.vertex_position_m(coordinate));
        }
    encoder.u64(topology.local_cell_count());
    for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
      encoder.u64(topology.global_cell_id(cell));
      append(encoder, geometry.cell_center_m(cell));
      encoder.f64(geometry.cell_volume_m3(cell));
      encoder.f64(geometry.minimum_jacobian_determinant_m3(cell));
      if (topology.cell_ownership(cell) == mesh::EntityOwnership::owned)
        append(encoder, geometry.cell_closure_m2(cell));
    }
    encoder.u64(topology.local_face_count());
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      encoder.u64(topology.global_face_id(face));
      append(encoder, geometry.face_center_m(face));
      append(encoder, geometry.face_displacement_m(face));
      append(encoder,
             geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
      append(encoder,
             geometry.face_area_vector_m2(face, mesh::FaceSide::neighbour));
      encoder.f64(geometry.face_area_m2(face));
      encoder.f64(geometry.face_skewness(face));
      encoder.f64(geometry.face_non_orthogonality_degrees(face));
    }
  });
}

std::uint64_t boundary_fingerprint(
    const boundary::BoundaryRegistry &boundaries) {
  return make_fingerprint("hundun.checkpoint-v2.boundary.v1",
                          [&](Encoder &encoder) {
    encoder.u64(boundaries.scalar_count());
    for (std::size_t scalar = 0; scalar < boundaries.scalar_count(); ++scalar)
      encoder.string(std::string(boundaries.scalar_name(scalar)));
    encoder.boolean(boundaries.open_domain());
    encoder.boolean(boundaries.velocity_inlet_patch_id().has_value());
    if (boundaries.velocity_inlet_patch_id())
      encoder.u32(*boundaries.velocity_inlet_patch_id());
    encoder.boolean(boundaries.pressure_outlet_patch_id().has_value());
    if (boundaries.pressure_outlet_patch_id())
      encoder.u32(*boundaries.pressure_outlet_patch_id());
    for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
      const auto &item = boundaries.patch(patch);
      encoder.u32(item.stable_id());
      encoder.string(std::string(item.name()));
      encoder.u8(code(item.kind()));
      encoder.u8(code(item.velocity_rule()));
      encoder.u8(code(item.pressure_rule()));
      encoder.u8(code(item.density_rule()));
      encoder.u8(code(item.enthalpy_rule()));
      encoder.u8(code(item.scalar_rule()));
      encoder.u8(code(item.mass_flux_rule()));
      encoder.boolean(item.paired_patch_id().has_value());
      if (item.paired_patch_id())
        encoder.u32(*item.paired_patch_id());
      encoder.boolean(item.inlet_state().has_value());
      if (item.inlet_state()) {
        append(encoder, item.inlet_state()->velocity_m_per_s);
        encoder.f64(item.inlet_state()->density_kg_per_m3);
        encoder.f64(item.inlet_state()->enthalpy_J_per_kg);
        append_optional(encoder, item.inlet_state()->temperature_K);
        encoder.u64(item.inlet_state()->scalar_values.size());
        for (const auto scalar : item.inlet_state()->scalar_values)
          encoder.f64(scalar);
      }
      append_optional(encoder, item.pressure_value_pa());
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
      encoder.u8(code(descriptor.space));
      encoder.u8(code(descriptor.scalar_type));
      encoder.u32(descriptor.components);
      encoder.i32(descriptor.ghost_width);
      encoder.boolean(descriptor.conservative);
      encoder.u8(code(descriptor.restart));
      encoder.u8(code(descriptor.output));
    }
  });
}

std::uint64_t local_layout_fingerprint(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, int rank, int rank_count) {
  return make_fingerprint("hundun.checkpoint-v2.local-layout.v1",
                          [&](Encoder &encoder) {
    const auto box = decomposition.owned_box();
    const auto extent = decomposition.local_extent();
    const auto grid = decomposition.process_grid();
    encoder.i32(rank);
    encoder.i32(rank_count);
    append(encoder, grid);
    append(encoder, box);
    append(encoder, extent);
    encoder.u64(topology.local_cell_count());
    for (std::size_t cell = 0; cell < topology.local_cell_count(); ++cell) {
      encoder.u8(topology.cell_ownership(cell) == mesh::EntityOwnership::owned
                     ? 0U : 1U);
      encoder.u64(topology.global_cell_id(cell));
    }
    encoder.u64(topology.local_face_count());
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      encoder.u8(topology.face_ownership(face) == mesh::EntityOwnership::owned
                     ? 0U : 1U);
      encoder.u64(topology.global_face_id(face));
    }
  });
}

std::array<std::uint64_t, 5> fingerprints(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::FlowCaseConfig &config,
    const runtime::FieldRegistry &registry, const FlowFieldIds &fields,
    std::uint64_t &collective_count) {
  const auto local_topology =
      local_topology_fingerprint(decomposition, topology, mpi.rank(),
                                 mpi.size());
  const auto local_geometry =
      local_geometry_fingerprint(topology, geometry, config, mpi.rank());
  const auto topology_parts = runtime::checkpoint_v2::allgather_u64(
      mpi, &local_topology, 1U, collective_count,
      "MPI_Allgather(Checkpoint topology fingerprints)");
  const auto geometry_parts = runtime::checkpoint_v2::allgather_u64(
      mpi, &local_geometry, 1U, collective_count,
      "MPI_Allgather(Checkpoint geometry fingerprints)");
  const auto topology_common =
      make_fingerprint("hundun.checkpoint-v2.topology-common.v1",
                       [&](Encoder &encoder) {
    encoder.i32(mpi.size());
    append(encoder, decomposition.process_grid());
    append(encoder, topology.global_extent());
    encoder.u64(topology_parts.size());
    for (int rank = 0; rank < mpi.size(); ++rank) {
      encoder.i32(rank);
      encoder.u64(topology_parts[static_cast<std::size_t>(rank)]);
    }
  });
  const auto geometry_common =
      make_fingerprint("hundun.checkpoint-v2.geometry-common.v1",
                       [&](Encoder &encoder) {
    encoder.u8(geometry.mapping_kind() == mesh::MappingKind::uniform_box
                   ? 0U : 1U);
    append(encoder, geometry.global_extent());
    append(encoder, geometry.origin_m());
    append(encoder, geometry.length_m());
    append_optional(encoder, geometry.uniform_spacing_m());
    append_optional(encoder, config.mesh.warp_amplitude);
    encoder.u64(geometry_parts.size());
    for (int rank = 0; rank < mpi.size(); ++rank) {
      encoder.i32(rank);
      encoder.u64(geometry_parts[static_cast<std::size_t>(rank)]);
    }
  });
  return {resolved_fingerprint(config), topology_common, geometry_common,
          boundary_fingerprint(boundaries),
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
  encoder.u8(code(metadata.order));
  encoder.u32(time.schema_version);
  encoder.u64(time.accepted_step);
  encoder.f64(time.proposed_next_dt_s);
  encoder.f64(time.last_accepted_dt_s);
  encoder.u8(code(time.last_accepted_order));
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
    encoder.u8(code(closure->mode));
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

std::uint64_t expected_global_payload_size(config::DensityModel model,
                                           bool open_domain) {
  if (model != config::DensityModel::ideal_gas)
    return 106U;
  return open_domain ? 124U : 132U;
}

GlobalPayload decode_global_payload(const ByteVector &bytes,
                                    std::uint64_t expected_size) {
  if (bytes.size() != runtime::checkpoint_v2::checked_size(expected_size))
    throw runtime::Error("Checkpoint v2 global payload size is invalid");
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
  result.metadata.order =
      order == 0U ? MomentumTimeOrder::backward_euler
                  : MomentumTimeOrder::bdf2;
  result.time.schema_version = decoder.u32();
  result.time.accepted_step = decoder.u64();
  result.time.proposed_next_dt_s = decoder.f64();
  result.time.last_accepted_dt_s = decoder.f64();
  const auto last_order = decoder.u8();
  if (last_order > 1U)
    throw runtime::Error("Checkpoint v2 controller order is invalid");
  result.time.last_accepted_order =
      last_order == 0U ? MomentumTimeOrder::backward_euler
                       : MomentumTimeOrder::bdf2;
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
    closure.mode = mode == 0U ? IdealGasPressureMode::closed_dynamic
                              : IdealGasPressureMode::open_fixed;
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

struct RankPayloadShape final {
  std::size_t cell_count{};
  std::size_t record_count{};
  std::uint64_t logical_bytes{};
  std::uint64_t payload_bytes{};
};

RankPayloadShape expected_rank_payload_shape(const FlowState &state) {
  const auto layout = detail::FlowStateCheckpointAccess::layout(state);
  if (layout.cell_interior_extent.x <= 0 ||
      layout.cell_interior_extent.y <= 0 ||
      layout.cell_interior_extent.z <= 0)
    throw runtime::Error("Checkpoint v2 cell extent is invalid");
  std::size_t cells = runtime::checkpoint_v2::checked_product(
      static_cast<std::size_t>(layout.cell_interior_extent.x),
      static_cast<std::size_t>(layout.cell_interior_extent.y));
  cells = runtime::checkpoint_v2::checked_product(
      cells, static_cast<std::size_t>(layout.cell_interior_extent.z));
  const auto transported = state.fields().transported_cell_fields.size();
  if (transported >
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max() / 2U -
                               5U))
    throw runtime::Error("Checkpoint v2 transported-field count is invalid");
  const auto records = runtime::checkpoint_v2::checked_product(
      2U, runtime::checkpoint_v2::checked_sum_u64(5U, transported));
  std::uint64_t one_layer_values{};
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values, static_cast<std::uint64_t>(cells));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values,
      runtime::checkpoint_v2::checked_product(cells, 3U));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values, static_cast<std::uint64_t>(cells));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values,
      runtime::checkpoint_v2::checked_product(layout.face_count, 3U));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values, static_cast<std::uint64_t>(layout.face_count));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values,
      runtime::checkpoint_v2::checked_product(cells, transported));
  const auto logical = runtime::checkpoint_v2::checked_product(
      runtime::checkpoint_v2::checked_product(one_layer_values, 2U),
      sizeof(double));
  const auto headers =
      runtime::checkpoint_v2::checked_product(records, 32U);
  const auto payload = runtime::checkpoint_v2::checked_sum_u64(
      32U, runtime::checkpoint_v2::checked_sum_u64(headers, logical));
  return {cells, records, logical, payload};
}

std::pair<FlowLayerValues, FlowLayerValues>
decode_rank_payload(const ByteVector &bytes, const FlowState &state,
                    std::uint64_t &logical_bytes) {
  const auto shape = expected_rank_payload_shape(state);
  if (bytes.size() !=
      runtime::checkpoint_v2::checked_size(shape.payload_bytes))
    throw runtime::Error("Checkpoint v2 rank payload size is invalid");
  Decoder decoder(bytes);
  const auto layout = detail::FlowStateCheckpointAccess::layout(state);
  const auto &fields = state.fields();
  if (decoder.u32() != 1U ||
      !same(runtime::Int3{decoder.i32(), decoder.i32(), decoder.i32()},
            layout.cell_interior_extent) ||
      decoder.u64() != layout.face_count ||
      decoder.u32() != fields.transported_cell_fields.size() ||
      decoder.u32() != shape.record_count)
    throw runtime::Error("Checkpoint v2 rank payload layout is invalid");
  FlowLayerValues layers[2];
  const std::size_t cell_count = shape.cell_count;
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
          byte_size != runtime::checkpoint_v2::checked_product(
                           runtime::checkpoint_v2::checked_size(items),
                           runtime::checkpoint_v2::checked_product(
                               components, sizeof(double))))
        throw runtime::Error("Checkpoint v2 field record is invalid");
      const auto value_count = runtime::checkpoint_v2::checked_product(
          runtime::checkpoint_v2::checked_size(items), components);
      std::vector<double> values;
      values.reserve(value_count);
      for (std::size_t index = 0; index < value_count; ++index)
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
  logical_bytes = shape.logical_bytes;
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

std::string rank_filename(int rank) {
  std::array<char, 32> buffer{};
  const int count = std::snprintf(buffer.data(), buffer.size(),
                                  "rank-%06d.v2.bin", rank);
  if (count <= 0 || static_cast<std::size_t>(count) >= buffer.size())
    throw runtime::Error("Checkpoint v2 rank filename is invalid");
  return std::string(buffer.data(), static_cast<std::size_t>(count));
}

using Convergence = runtime::checkpoint_v2::CollectiveResult;
Convergence path_agrees(const runtime::MpiContext &mpi,
                        const std::filesystem::path &path,
                        std::uint64_t &collective_count) {
  const auto text = path.lexically_normal().generic_string();
  bool local_valid = !text.empty() && text.find('\0') == std::string::npos;
  Encoder encoder;
  try {
    if (local_valid)
      encoder.string(text);
  } catch (...) {
    local_valid = false;
  }
  const auto agreement = runtime::checkpoint_v2::opaque_bytes_agreement(
      mpi, encoder.bytes(), collective_count,
      "MPI_Allreduce(Checkpoint v2 path agreement)");
  return runtime::checkpoint_v2::converge_phase(
      mpi, local_valid && agreement.ok, collective_count,
      "MPI_Allreduce(Checkpoint v2 path validity)");
}

Convergence converge(const runtime::MpiContext &mpi, bool local_ok,
                     std::uint64_t &collective_count,
                     std::string_view operation) {
  return runtime::checkpoint_v2::converge_phase(
      mpi, local_ok, collective_count, operation);
}

std::uint64_t common_fingerprint(
    const std::array<std::uint64_t, 5> &items, int ranks, runtime::Int3 grid,
    const ByteVector &global) {
  Encoder encoder;
    for (auto item : items)
      encoder.u64(item);
    encoder.i32(ranks);
    encoder.i32(grid.x);
    encoder.i32(grid.y);
    encoder.i32(grid.z);
    encoder.raw(global.data(), global.size());
  return runtime::checkpoint_v2::crc64_ecma(encoder.bytes().data(),
                                             encoder.bytes().size());
}

bool directory_inventory_valid(const std::filesystem::path &directory,
                               int ranks) {
  std::vector<std::string> expected{"COMPLETED", "manifest.v2.bin"};
  for (int rank = 0; rank < ranks; ++rank)
    expected.push_back(rank_filename(rank));
  return runtime::checkpoint_v2::exact_directory_inventory(directory,
                                                            expected);
}

CheckpointV2FailureReason file_failure_reason(
    runtime::checkpoint_v2::NumericFileFailure failure) noexcept {
  return failure == runtime::checkpoint_v2::NumericFileFailure::filesystem
             ? CheckpointV2FailureReason::filesystem
             : CheckpointV2FailureReason::file_integrity;
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
CheckpointV2ReadResult::CheckpointV2ReadResult(
    CheckpointV2Report report) noexcept
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
    std::optional<IdealGasClosureState> closure, bool restored) noexcept {
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
    std::uint64_t time_bits{};
    std::memcpy(&time_bits, &values.time_s, sizeof(time_bits));
    const std::array<std::uint64_t, 3> local_common{
        static_cast<std::uint8_t>(reason), values.step, time_bits};
    const auto common = runtime::checkpoint_v2::allgather_u64(
        mpi, local_common.data(), local_common.size(), collectives,
        "MPI_Allgather(Checkpoint failure report)");
    const auto failing_offset =
        static_cast<std::size_t>(rank) * local_common.size();
    reason = static_cast<CheckpointV2FailureReason>(
        common[failing_offset]);
    values.step = common[1U];
    std::memcpy(&values.time_s, &common[2U], sizeof(values.time_s));
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
  const auto path_status = path_agrees(mpi, directory, collectives);
  if (!path_status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, path_status.failing_rank);
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
    if (!local_ok)
      throw runtime::Error("Checkpoint v2 preflight identity is invalid");
    metadata = state.metadata();
    values.step = metadata.step;
    values.time_s = metadata.time_s;
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
    identity = fingerprints(mpi, decomposition, topology, geometry, boundaries,
                            config, registry, state.fields(), collectives);
    local_layout =
        local_layout_fingerprint(decomposition, topology, mpi.rank(),
                                 mpi.size());
    global_payload = encode_global_payload(metadata, time, closure);
    if (global_payload.size() !=
        expected_global_payload_size(config.density_model,
                                     boundaries.open_domain()))
      throw runtime::Error("Checkpoint v2 global payload size is invalid");
    rank_payload = encode_rank_payload(state, values.local_logical_bytes);
    const auto rank_shape = expected_rank_payload_shape(state);
    if (rank_payload.size() != rank_shape.payload_bytes ||
        values.local_logical_bytes != rank_shape.logical_bytes)
      throw runtime::Error("Checkpoint v2 rank payload size is invalid");
    rank_wrapper = runtime::checkpoint_v2::encode_rank_wrapper(
        mpi.rank(), mpi.size(), rank_payload);
    values.local_actual_bytes = rank_wrapper.size();
    values.local_crc64 = runtime::checkpoint_v2::crc64_ecma(
        rank_wrapper.data(), rank_wrapper.size());
    owned = decomposition.owned_box();
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    local_ok = false;
  }
  auto status =
      converge(mpi, local_ok, collectives, "MPI_Allreduce(Checkpoint preflight)");
  if (!status.ok) {
    values.fingerprint = CheckpointV2CheckStatus::not_checked;
    values.partition = local_reason == CheckpointV2FailureReason::layout
                           ? CheckpointV2CheckStatus::failed
                           : CheckpointV2CheckStatus::not_checked;
    return fail(local_reason, CheckpointV2Phase::preflight,
                status.failing_rank);
  }
  Encoder common_authority;
  for (const auto fingerprint : identity)
    common_authority.u64(fingerprint);
  common_authority.i32(mpi.size());
  append(common_authority, decomposition.process_grid());
  common_authority.u64(global_payload.size());
  common_authority.raw(global_payload.data(), global_payload.size());
  status = runtime::checkpoint_v2::opaque_bytes_agreement(
      mpi, common_authority.bytes(), collectives,
      "MPI_Allreduce(Checkpoint common authority)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, status.failing_rank);
  values.step = metadata.step;
  values.time_s = metadata.time_s;
  values.fingerprint = CheckpointV2CheckStatus::passed;
  values.partition = CheckpointV2CheckStatus::passed;

  if (config.density_model == config::DensityModel::ideal_gas) {
    const bool model_valid = detail::validate_ideal_gas_restore_state(
        mpi, topology, geometry, boundaries,
        config.physics.cp_J_per_kg_K.value_or(
            std::numeric_limits<double>::quiet_NaN()),
        config.physics.gas_constant_J_per_kg_K.value_or(
            std::numeric_limits<double>::quiet_NaN()),
        config.physics.thermodynamic_pressure_pa.value_or(
            std::numeric_limits<double>::quiet_NaN()),
        history, committed, *closure, collectives);
    status = converge(mpi, model_valid, collectives,
                      "MPI_Allreduce(Checkpoint ideal-gas write state)");
    if (!status.ok)
      return fail(CheckpointV2FailureReason::state,
                  CheckpointV2Phase::rank_payload, status.failing_rank);
  }

  bool filesystem_ok = true;
  if (mpi.rank() == 0) {
    try {
      runtime::checkpoint_v2::create_directory_exclusive(directory);
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
    values.crc_check_count = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, rank_written ? 1U : 0U, collectives,
        "MPI_Allreduce(Checkpoint verified rank files)");
    values.rank_crc = rank_written ? CheckpointV2CheckStatus::passed
                                   : CheckpointV2CheckStatus::not_checked;
    values.exact_size_eof = CheckpointV2CheckStatus::not_checked;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::rank_temporary_file, status.failing_rank);
  }
  values.rank_crc = CheckpointV2CheckStatus::passed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size());

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
    values.file_count = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, rank_published ? 1U : 0U, collectives,
        "MPI_Allreduce(Checkpoint published rank files)");
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::rank_publish, status.failing_rank);
  }
  values.file_count = static_cast<std::uint64_t>(mpi.size());

  const std::array<std::uint64_t, 11> local_record{
      static_cast<std::uint32_t>(mpi.rank()),
      static_cast<std::uint32_t>(owned.begin.x),
      static_cast<std::uint32_t>(owned.begin.y),
      static_cast<std::uint32_t>(owned.begin.z),
      static_cast<std::uint32_t>(owned.end.x),
      static_cast<std::uint32_t>(owned.end.y),
      static_cast<std::uint32_t>(owned.end.z),
      values.local_logical_bytes, values.local_actual_bytes,
      values.local_crc64, local_layout};
  const auto gathered_records = runtime::checkpoint_v2::allgather_u64(
      mpi, local_record.data(), local_record.size(), collectives,
      "MPI_Allgather(Checkpoint rank records)");
  const auto global_logical = runtime::checkpoint_v2::allreduce_sum_u64(
      mpi, values.local_logical_bytes, collectives,
      "MPI_Allreduce(Checkpoint logical bytes)");
  const auto rank_actual = runtime::checkpoint_v2::allreduce_sum_u64(
      mpi, values.local_actual_bytes, collectives,
      "MPI_Allreduce(Checkpoint actual bytes)");

  std::uint64_t manifest_size = 0U;
  std::uint64_t manifest_crc = 0U;
  std::uint64_t marker_size = 0U;
  bool manifest_ok = true;
  bool manifest_verified = false;
  if (mpi.rank() == 0) {
    try {
      runtime::checkpoint_v2::Manifest manifest;
      manifest.rank_count = static_cast<std::uint32_t>(mpi.size());
      manifest.process_grid = decomposition.process_grid();
      manifest.fingerprints = identity;
      manifest.global_payload = global_payload;
      for (int rank = 0; rank < mpi.size(); ++rank) {
        const auto offset =
            static_cast<std::size_t>(rank) * local_record.size();
        const auto signed32 = [&](std::size_t index) {
          const auto bits =
              static_cast<std::uint32_t>(gathered_records[offset + index]);
          std::int32_t result{};
          std::memcpy(&result, &bits, sizeof(result));
          return result;
        };
        manifest.ranks.push_back(
            {signed32(0U),
             {signed32(1U), signed32(2U), signed32(3U)},
             {signed32(4U), signed32(5U), signed32(6U)},
             rank_filename(rank), gathered_records[offset + 7U],
             gathered_records[offset + 8U], gathered_records[offset + 9U],
             gathered_records[offset + 10U]});
      }
      const auto bytes = runtime::checkpoint_v2::encode_manifest(manifest);
      const auto expected_size = runtime::checkpoint_v2::checked_sum_u64(
          runtime::checkpoint_v2::checked_sum_u64(
              84U, static_cast<std::uint64_t>(global_payload.size())),
          runtime::checkpoint_v2::checked_product(
              static_cast<std::size_t>(mpi.size()), 82U));
      if (bytes.size() != expected_size)
        throw runtime::Error("Checkpoint v2 manifest size is invalid");
      const auto temp = directory / "manifest.v2.bin.tmp";
      const auto verified =
          runtime::checkpoint_v2::write_verified_temporary(temp, bytes);
      manifest_size = verified.actual_size;
      manifest_crc = verified.crc64;
      manifest_verified = true;
      runtime::checkpoint_v2::publish_no_overwrite(
          temp, directory / "manifest.v2.bin");
    } catch (...) {
      manifest_ok = false;
    }
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest)");
  const bool manifest_verified_common =
      runtime::checkpoint_v2::allreduce_sum_u64(
          mpi, manifest_verified ? 1U : 0U, collectives,
          "MPI_Allreduce(Checkpoint manifest verification)") == 1U;
  if (!status.ok) {
    values.manifest_crc = manifest_verified_common
                              ? CheckpointV2CheckStatus::passed
                              : CheckpointV2CheckStatus::not_checked;
    if (manifest_verified_common)
      ++values.crc_check_count;
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::manifest, status.failing_rank);
  }
  const std::array<std::uint64_t, 2> local_manifest{
      mpi.rank() == 0 ? manifest_size : 0U,
      mpi.rank() == 0 ? manifest_crc : 0U};
  const auto manifest_authority = runtime::checkpoint_v2::allgather_u64(
      mpi, local_manifest.data(), local_manifest.size(), collectives,
      "MPI_Allgather(Checkpoint manifest authority)");
  manifest_size = manifest_authority[0U];
  manifest_crc = manifest_authority[1U];
  values.manifest_crc64 = manifest_crc;
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 1U;

  bool marker_ok = true;
  bool marker_verified = false;
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
      marker_verified = true;
      runtime::checkpoint_v2::publish_no_overwrite(temp,
                                                    directory / "COMPLETED");
    } catch (...) {
      marker_ok = false;
    }
  }
  status = converge(mpi, marker_ok, collectives,
                    "MPI_Allreduce(Checkpoint marker)");
  const bool marker_verified_common =
      runtime::checkpoint_v2::allreduce_sum_u64(
          mpi, marker_verified ? 1U : 0U, collectives,
          "MPI_Allreduce(Checkpoint marker verification)") == 1U;
  if (!status.ok) {
    if (marker_verified_common)
      ++values.crc_check_count;
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::filesystem,
                CheckpointV2Phase::completed_marker, status.failing_rank);
  }
  const auto marker_authority = runtime::checkpoint_v2::allgather_u64(
      mpi, &marker_size, 1U, collectives,
      "MPI_Allgather(Checkpoint marker size)");
  marker_size = marker_authority[0U];
  ++values.crc_check_count;
  values.global_logical_bytes = global_logical;
  values.global_actual_bytes = runtime::checkpoint_v2::checked_sum_u64(
      runtime::checkpoint_v2::checked_sum_u64(rank_actual, manifest_size),
      marker_size);
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
    std::uint64_t time_bits{};
    std::memcpy(&time_bits, &values.time_s, sizeof(time_bits));
    const std::array<std::uint64_t, 3> local_common{
        static_cast<std::uint8_t>(reason), values.step, time_bits};
    const auto common = runtime::checkpoint_v2::allgather_u64(
        mpi, local_common.data(), local_common.size(), collectives,
        "MPI_Allgather(Checkpoint read failure report)");
    const auto failing_offset =
        static_cast<std::size_t>(rank) * local_common.size();
    reason = static_cast<CheckpointV2FailureReason>(
        common[failing_offset]);
    values.step = common[1U];
    std::memcpy(&values.time_s, &common[2U], sizeof(values.time_s));
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
  const auto path_status = path_agrees(mpi, directory, collectives);
  if (!path_status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, path_status.failing_rank);
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
  bool marker_exact = false;
  bool marker_exact_failed = false;
  auto marker_reason = CheckpointV2FailureReason::file_integrity;
  try {
    const auto marker_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "COMPLETED", 40U);
    marker_exact = true;
    marker =
        runtime::checkpoint_v2::decode_completed_marker(marker_bytes);
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    marker_reason = file_failure_reason(error.failure());
    marker_exact_failed =
        error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    marker_ok = false;
  } catch (...) {
    marker_ok = false;
  }
  status = converge(mpi, marker_ok, collectives,
                    "MPI_Allreduce(Checkpoint marker read)");
  if (!status.ok) {
    values.exact_size_eof =
        marker_exact ? CheckpointV2CheckStatus::passed
                     : marker_exact_failed
                           ? CheckpointV2CheckStatus::failed
                           : CheckpointV2CheckStatus::not_checked;
    return fail(marker_reason,
                CheckpointV2Phase::marker_read, status.failing_rank);
  }
  values.file_count = 1U;

  bool manifest_ok = true;
  bool manifest_exact = false;
  bool manifest_exact_failed = false;
  bool manifest_crc_checked = false;
  bool manifest_crc_match = false;
  auto manifest_reason = CheckpointV2FailureReason::file_integrity;
  try {
    const auto global_size = expected_global_payload_size(
        config.density_model, boundaries.open_domain());
    const auto expected_manifest_size =
        runtime::checkpoint_v2::checked_sum_u64(
            runtime::checkpoint_v2::checked_sum_u64(84U, global_size),
            runtime::checkpoint_v2::checked_product(
                static_cast<std::size_t>(mpi.size()), 82U));
    if (marker.manifest_actual_size != expected_manifest_size)
      manifest_exact_failed = true;
    if (manifest_exact_failed)
      throw runtime::Error("Checkpoint v2 manifest declared size is invalid");
    manifest_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "manifest.v2.bin", marker.manifest_actual_size);
    manifest_exact = true;
    values.manifest_crc64 = runtime::checkpoint_v2::crc64_ecma(
        manifest_bytes.data(), manifest_bytes.size());
    manifest_crc_checked = true;
    if (values.manifest_crc64 != marker.manifest_crc64)
      throw runtime::Error("Checkpoint v2 manifest CRC is invalid");
    manifest_crc_match = true;
    manifest = runtime::checkpoint_v2::decode_manifest(
        manifest_bytes, static_cast<std::uint32_t>(mpi.size()), global_size);
    global = decode_global_payload(manifest.global_payload, global_size);
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    manifest_reason = file_failure_reason(error.failure());
    manifest_exact_failed =
        error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    manifest_ok = false;
  } catch (...) {
    manifest_ok = false;
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest read)");
  if (!status.ok) {
    values.manifest_crc =
        !manifest_crc_checked ? CheckpointV2CheckStatus::not_checked
                              : manifest_crc_match
                                    ? CheckpointV2CheckStatus::passed
                                    : CheckpointV2CheckStatus::failed;
    values.exact_size_eof =
        manifest_exact ? CheckpointV2CheckStatus::passed
                       : manifest_exact_failed
                             ? CheckpointV2CheckStatus::failed
                             : CheckpointV2CheckStatus::not_checked;
    values.crc_check_count = manifest_crc_match ? 1U : 0U;
    return fail(manifest_reason,
                CheckpointV2Phase::manifest_read, status.failing_rank);
  }
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  values.file_count = 2U;
  values.crc_check_count = 1U;
  values.step = global.metadata.step;
  values.time_s = global.metadata.time_s;

  bool inventory_ok = false;
  auto inventory_reason = CheckpointV2FailureReason::file_integrity;
  try {
    inventory_ok = directory_inventory_valid(directory, mpi.size());
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    inventory_reason = file_failure_reason(error.failure());
    inventory_ok = false;
  } catch (...) {
    inventory_ok = false;
  }
  status = converge(mpi, inventory_ok, collectives,
                    "MPI_Allreduce(Checkpoint inventory)");
  if (!status.ok)
    return fail(inventory_reason,
                CheckpointV2Phase::manifest_read, status.failing_rank);

  bool partition_ok = true;
  try {
    const auto &record =
        manifest.ranks.at(static_cast<std::size_t>(mpi.rank()));
    partition_ok =
        manifest.rank_count == static_cast<std::uint32_t>(mpi.size()) &&
        same(manifest.process_grid, decomposition.process_grid()) &&
        record.rank == mpi.rank() &&
        same(runtime::Box3{record.owned_box_begin, record.owned_box_end},
             decomposition.owned_box()) &&
        record.local_layout_fingerprint ==
            local_layout_fingerprint(decomposition, topology, mpi.rank(),
                                     mpi.size());
  } catch (...) {
    partition_ok = false;
  }
  status = converge(mpi, partition_ok, collectives,
                    "MPI_Allreduce(Checkpoint partition)");
  if (!status.ok) {
    values.partition = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::layout,
                CheckpointV2Phase::manifest_read, status.failing_rank);
  }
  values.partition = CheckpointV2CheckStatus::passed;

  bool fingerprint_ok = true;
  try {
    const auto &registry = detail::FlowStateCheckpointAccess::registry(state);
    const auto expected = fingerprints(
        mpi, decomposition, topology, geometry, boundaries, config, registry,
        state.fields(), collectives);
    fingerprint_ok =
        manifest.fingerprints == expected &&
        marker.common_fingerprint ==
            common_fingerprint(expected, mpi.size(),
                               decomposition.process_grid(),
                               manifest.global_payload);
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    fingerprint_ok = false;
  }
  status = converge(mpi, fingerprint_ok, collectives,
                    "MPI_Allreduce(Checkpoint fingerprints)");
  if (!status.ok) {
    values.fingerprint = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::file_integrity,
                CheckpointV2Phase::manifest_read, status.failing_rank);
  }
  values.fingerprint = CheckpointV2CheckStatus::passed;

  const bool global_state_ok =
      detail::TimeControlStateCodec::semantically_valid(
          config.time, config.density_model, global.metadata, global.time) &&
      (global.closure.has_value() ==
       (config.density_model == config::DensityModel::ideal_gas));
  status = converge(mpi, global_state_ok, collectives,
                    "MPI_Allreduce(Checkpoint global state)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare, status.failing_rank);

  ByteVector rank_bytes;
  std::pair<FlowLayerValues, FlowLayerValues> layers;
  bool rank_ok = true;
  bool rank_exact = false;
  bool rank_exact_failed = false;
  bool rank_crc_checked = false;
  bool rank_crc_match = false;
  auto rank_reason = CheckpointV2FailureReason::file_integrity;
  try {
    const auto shape = expected_rank_payload_shape(state);
    const auto &record =
        manifest.ranks.at(static_cast<std::size_t>(mpi.rank()));
    const auto expected_actual =
        runtime::checkpoint_v2::checked_sum_u64(32U, shape.payload_bytes);
    if (record.actual_byte_size != expected_actual ||
        record.logical_byte_size != shape.logical_bytes) {
      rank_exact_failed = record.actual_byte_size != expected_actual;
      throw runtime::Error("Checkpoint v2 rank record sizes are invalid");
    }
    rank_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / record.filename, expected_actual);
    rank_exact = true;
    values.local_actual_bytes = rank_bytes.size();
    values.local_crc64 = runtime::checkpoint_v2::crc64_ecma(
        rank_bytes.data(), rank_bytes.size());
    rank_crc_checked = true;
    if (values.local_crc64 != record.crc64)
      throw runtime::Error("Checkpoint v2 rank CRC is invalid");
    rank_crc_match = true;
    const auto wrapper =
        runtime::checkpoint_v2::decode_rank_wrapper(rank_bytes,
                                                    shape.payload_bytes);
    if (wrapper.rank != mpi.rank() || wrapper.rank_count != mpi.size())
      throw runtime::Error("Checkpoint v2 rank wrapper identity is invalid");
    std::uint64_t decoded_logical{};
    layers = decode_rank_payload(wrapper.payload, state, decoded_logical);
    if (decoded_logical != record.logical_byte_size)
      throw runtime::Error("Checkpoint v2 rank logical size is invalid");
    values.local_logical_bytes = decoded_logical;
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    rank_reason = file_failure_reason(error.failure());
    rank_exact_failed =
        error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    rank_ok = false;
  } catch (...) {
    rank_ok = false;
  }
  status = converge(mpi, rank_ok, collectives,
                    "MPI_Allreduce(Checkpoint rank read)");
  if (!status.ok) {
    const auto exact_count = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, rank_exact ? 1U : 0U, collectives,
        "MPI_Allreduce(Checkpoint exact rank files)");
    const auto exact_failure_count =
        runtime::checkpoint_v2::allreduce_sum_u64(
            mpi, rank_exact_failed ? 1U : 0U, collectives,
            "MPI_Allreduce(Checkpoint invalid rank sizes)");
    const auto crc_count = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, rank_crc_match ? 1U : 0U, collectives,
        "MPI_Allreduce(Checkpoint matching rank CRCs)");
    const auto complete_count = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, rank_ok ? 1U : 0U, collectives,
        "MPI_Allreduce(Checkpoint complete rank files)");
    values.file_count = runtime::checkpoint_v2::checked_sum_u64(
        2U, complete_count);
    values.crc_check_count = runtime::checkpoint_v2::checked_sum_u64(
        1U, crc_count);
    values.rank_crc =
        !rank_crc_checked ? CheckpointV2CheckStatus::not_checked
                          : rank_crc_match
                                ? CheckpointV2CheckStatus::passed
                                : CheckpointV2CheckStatus::failed;
    values.exact_size_eof =
        exact_count == static_cast<std::uint64_t>(mpi.size())
            ? CheckpointV2CheckStatus::passed
            : exact_failure_count != 0U
                  ? CheckpointV2CheckStatus::failed
                  : CheckpointV2CheckStatus::not_checked;
    return fail(rank_reason,
                CheckpointV2Phase::rank_read, status.failing_rank);
  }
  values.rank_crc = CheckpointV2CheckStatus::passed;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;

  bool state_ok =
      valid_layer(layers.first, config) &&
      valid_layer(layers.second, config);
  if (config.density_model == config::DensityModel::ideal_gas &&
      global.closure)
    state_ok =
        detail::validate_ideal_gas_restore_state(
            mpi, topology, geometry, boundaries,
            config.physics.cp_J_per_kg_K.value_or(
                std::numeric_limits<double>::quiet_NaN()),
            config.physics.gas_constant_J_per_kg_K.value_or(
                std::numeric_limits<double>::quiet_NaN()),
            config.physics.thermodynamic_pressure_pa.value_or(
                std::numeric_limits<double>::quiet_NaN()),
            layers.first, layers.second, *global.closure, collectives) &&
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

  values.publication = CheckpointV2CheckStatus::passed;
  values.rollback = CheckpointV2CheckStatus::not_checked;
  values.disposition = CheckpointV2Disposition::completed;
  values.reason = CheckpointV2FailureReason::none;
  values.phase = CheckpointV2Phase::restore_publish;
  values.lowest_failing_rank = -1;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
  values.global_logical_bytes = runtime::checkpoint_v2::allreduce_sum_u64(
      mpi, values.local_logical_bytes, collectives,
      "MPI_Allreduce(Checkpoint restored logical bytes)");
  const auto rank_actual = runtime::checkpoint_v2::allreduce_sum_u64(
      mpi, values.local_actual_bytes, collectives,
      "MPI_Allreduce(Checkpoint restored actual bytes)");
  values.global_actual_bytes = runtime::checkpoint_v2::checked_sum_u64(
      rank_actual, runtime::checkpoint_v2::checked_sum_u64(
                       marker.manifest_actual_size, 40U));
  values.collective_count =
      runtime::checkpoint_v2::checked_sum_u64(collectives, 1U);
  auto completed_report = detail::CheckpointV2Access::make(values);
  status = converge(mpi, true, collectives,
                    "MPI_Allreduce(Checkpoint restore success boundary)");
  if (!status.ok)
    throw runtime::Error("Checkpoint v2 success boundary is inconsistent");
  detail::FlowStateCheckpointAccess::publish_replacement(
      state, std::move(*replacement));
  return detail::CheckpointV2Access::make_read(
      std::move(completed_report), global.time, global.closure, true);
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

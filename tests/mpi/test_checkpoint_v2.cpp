// SPDX-License-Identifier: Apache-2.0

#include "checkpoint_v2_protocol.hpp"
#include "checkpoint_v2_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
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
#include "tests/support/test_main.hpp"

#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Bytes = std::vector<std::uint8_t>;

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
Bytes canonical_prefix(std::string_view domain) {
  Bytes bytes(domain.begin(), domain.end());
  append_u8(bytes, 0U);
  append_u32(bytes, 1U);
  return bytes;
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
    append_real3(bytes, owner);
    append_real3(bytes, {-owner.x, -owner.y, -owner.z});
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

void require_common_report(const hundun::runtime::MpiContext &mpi,
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

void run(const hundun::runtime::MpiContext &mpi) {
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
  auto make_state = [&] {
    return hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields, metadata);
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
  const auto independent_mesh_fingerprints =
      independent_common_mesh_fingerprints(mpi, decomposition, topology,
                                           geometry, config);
  if (mpi.rank() == 0) {
    std::ifstream manifest_stream(directory / "manifest.v2.bin",
                                  std::ios::binary);
    const Bytes manifest_bytes{std::istreambuf_iterator<char>(manifest_stream),
                               std::istreambuf_iterator<char>()};
    HUNDUN_CHECK(manifest_bytes.size() >= 72U);
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 40U) ==
                 independent_mesh_fingerprints.first);
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 48U) ==
                 independent_mesh_fingerprints.second);
    HUNDUN_CHECK(read_u64_at(manifest_bytes, 64U) ==
                 independent_field_schema_fingerprint(registry, fields));
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
  for (const auto &mutate : descriptor_mutations) {
    std::array<hundun::runtime::FieldDescriptor, 7> descriptors;
    for (std::size_t index = 0; index < descriptors.size(); ++index)
      descriptors[index] = registry.descriptor(ordered_schema_ids[index]);
    mutate(descriptors[0]);
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
    HUNDUN_CHECK(independent_field_schema_fingerprint(changed_registry,
                                                      changed_fields) !=
                 baseline_schema_fingerprint);
  }

  auto destination = make_state();
  auto different = history;
  for (double &item : different.density)
    item += 9.0;
  destination.seed_accepted_layers(different, different);
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

  using CheckpointAccess = hundun::flow::test::CheckpointV2TestAccess;
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
  for (std::size_t mutation_index = 0U;
       mutation_index < resolved_tuple_mutations.size(); ++mutation_index) {
    auto changed = config;
    resolved_tuple_mutations[mutation_index](changed);
    auto fingerprint_destination = make_state();
    fingerprint_destination.seed_accepted_layers(different, different);
    const auto before = CheckpointAccess::snapshot(fingerprint_destination);
    const auto rejected = hundun::flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, changed,
        fingerprint_destination, directory);
    if (rejected.restored())
      throw std::runtime_error(
          "Task23 resolved tuple mutation was not rejected: " +
          std::to_string(mutation_index));
    const auto after = CheckpointAccess::snapshot(fingerprint_destination);
    HUNDUN_CHECK(before.storage_bits == after.storage_bits);
    HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
        before.metadata, after.metadata));
  }

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
        fields, metadata);
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
    HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
            partition_before, CheckpointAccess::snapshot(partition_state)));
    for (const auto &view : partition_views)
      HUNDUN_CHECK(rejects([&] { static_cast<void>(view(0, 0, 0, 0)); }));
  }

  hundun::mesh::MeshGeometry changed_geometry(
      topology, hundun::mesh::AnalyticWarpedBoxMapping(config.mesh.origin_m,
                                                       config.mesh.length_m,
                                                       {0.01, -0.005, 0.0025}));
  auto geometry_destination = make_state();
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
  ideal_config.scalars.clear();
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
  ideal_fields.transported_cell_fields = {rho_h};
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
      std::vector<double>(cells, 101325.0 * 1000.0 / 287.05)};
  ideal_state.seed_accepted_layers(ideal_values, ideal_values);
  const hundun::flow::IdealGasClosureSpec ideal_spec{rho_h, 1000.0, 287.05,
                                                     101325.0};
  auto initial_closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, ideal_boundaries, mpi, ideal_registry, ideal_fields,
      ideal_state, ideal_spec);
  auto persisted_closure = initial_closure.state();
  persisted_closure.revision = 4U;
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
  if (mpi.rank() == 0)
    std::filesystem::remove_all(ideal_directory);
  mpi.barrier();
  const auto ideal_written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, ideal_boundaries, ideal_config,
      ideal_state, ideal_controller_state, persisted_closure, ideal_directory);
  HUNDUN_CHECK(ideal_written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
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
      std::vector<hundun::config::InletScalarValue>{};
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
  open_fields.transported_cell_fields = {open_rho_h};
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
  open_values.transported_cell_fields = {std::vector<double>(
      open_topology.owned_cell_count(), ideal_density * 300000.0)};
  const auto make_open_state = [&] {
    auto result = hundun::flow::FlowState::create(
        open_registry,
        {open_decomposition.local_extent(), open_topology.local_face_count()},
        open_fields, metadata);
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
    require_common_report(mpi, failed_phase);
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
  const auto continuation_directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-continuation-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(continuation_directory);
  mpi.barrier();
  const auto continuation_write = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, constant_boundaries,
      constant_config, uninterrupted, uninterrupted_controller.state(),
      std::nullopt, continuation_directory);
  HUNDUN_CHECK(continuation_write.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);

  struct ReadPhaseCase final {
    std::uint32_t calls_before;
    hundun::flow::CheckpointV2Phase phase;
  };
  constexpr std::array read_phase_cases{
      ReadPhaseCase{0U, hundun::flow::CheckpointV2Phase::marker_read},
      ReadPhaseCase{1U, hundun::flow::CheckpointV2Phase::manifest_read},
      ReadPhaseCase{2U, hundun::flow::CheckpointV2Phase::rank_read}};
  for (const auto &phase_case : read_phase_cases) {
    auto failed_destination = make_state();
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
    require_common_report(mpi, failed_phase.report());
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
    require_common_report(mpi, invalid_read.report());
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
  auto invalid_generic_destination = make_state();
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
  HUNDUN_CHECK(hundun::flow::test::checkpoint_v2_failed_read_preserved_values(
      before_invalid_generic,
      CheckpointAccess::snapshot(invalid_generic_destination)));
  require_common_report(mpi, invalid_generic_read.report());
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
    require_common_report(mpi, invalid_ideal_read.report());
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
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi);
  });
}

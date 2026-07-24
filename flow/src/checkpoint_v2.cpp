// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/checkpoint_v2.hpp"
#include "adaptive_time_control_detail.hpp"
#include "checkpoint_v2_detail.hpp"
#include "checkpoint_v2_protocol.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "checkpoint_v2_test_access.hpp"
#endif

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
#include <new>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

enum class CheckpointPreparationPoint : std::uint8_t {
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
  final_success_boundary
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
CheckpointPreparationPoint checkpoint_preparation_fault{
    CheckpointPreparationPoint::none};
std::uint32_t checkpoint_preparation_fault_calls_before{};

void inject_checkpoint_preparation_fault(CheckpointPreparationPoint point) {
  if (checkpoint_preparation_fault != point)
    return;
  if (checkpoint_preparation_fault_calls_before != 0U) {
    --checkpoint_preparation_fault_calls_before;
    return;
  }
  checkpoint_preparation_fault = CheckpointPreparationPoint::none;
  throw std::bad_alloc();
}
#else
void inject_checkpoint_preparation_fault(CheckpointPreparationPoint) {}
#endif

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

constexpr std::string_view kReportSealDomain =
    "hundun.checkpoint-v2.report-semantic.v1";
constexpr std::size_t kReportSealEncodedSize =
    sizeof(std::uint32_t) + kReportSealDomain.size() +
    4U * sizeof(std::uint8_t) + 2U * sizeof(std::int32_t) +
    11U * sizeof(std::uint64_t) + 8U * sizeof(std::uint8_t);
static_assert(kReportSealEncodedSize == 151U);

class ReportSealEncoder final {
public:
  void u8(std::uint8_t item) noexcept {
    if (size_ >= bytes_.size()) {
      valid_ = false;
      return;
    }
    bytes_[size_++] = item;
  }
  void u32(std::uint32_t item) noexcept {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      u8(static_cast<std::uint8_t>(item >> shift));
  }
  void i32(std::int32_t item) noexcept {
    std::uint32_t wire{};
    std::memcpy(&wire, &item, sizeof(wire));
    u32(wire);
  }
  void u64(std::uint64_t item) noexcept {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      u8(static_cast<std::uint8_t>(item >> shift));
  }
  void f64(double item) noexcept {
    std::uint64_t wire{};
    std::memcpy(&wire, &item, sizeof(wire));
    u64(wire);
  }
  void string(std::string_view item) noexcept {
    if (item.size() > std::numeric_limits<std::uint32_t>::max()) {
      valid_ = false;
      return;
    }
    u32(static_cast<std::uint32_t>(item.size()));
    if (!valid_ || item.size() > bytes_.size() - size_) {
      valid_ = false;
      return;
    }
    std::memcpy(bytes_.data() + size_, item.data(), item.size());
    size_ += item.size();
  }
  const std::uint8_t *data() const noexcept { return bytes_.data(); }
  std::size_t size() const noexcept { return size_; }
  bool valid() const noexcept { return valid_; }

private:
  std::array<std::uint8_t, kReportSealEncodedSize> bytes_{};
  std::size_t size_{};
  bool valid_{true};
};

std::uint64_t report_seal(const CheckpointV2Report &report) noexcept {
  ReportSealEncoder fixed;
  fixed.string(kReportSealDomain);
  fixed.u8(value(report.operation()));
  fixed.u8(value(report.disposition()));
  fixed.u8(value(report.reason()));
  fixed.u8(value(report.phase()));
  fixed.i32(report.rank());
  fixed.i32(report.lowest_failing_rank());
  fixed.u64(report.step());
  fixed.f64(report.time_s());
  fixed.u64(report.local_logical_bytes());
  fixed.u64(report.local_actual_bytes());
  fixed.u64(report.global_logical_bytes());
  fixed.u64(report.global_actual_bytes());
  fixed.u64(report.local_crc64());
  fixed.u64(report.manifest_crc64());
  fixed.u64(report.file_count());
  fixed.u64(report.crc_check_count());
  fixed.u64(report.collective_count());
  fixed.u8(value(report.rank_crc_status()));
  fixed.u8(value(report.manifest_crc_status()));
  fixed.u8(value(report.exact_size_and_eof_status()));
  fixed.u8(value(report.fingerprint_status()));
  fixed.u8(value(report.partition_status()));
  fixed.u8(value(report.transaction_entry_status()));
  fixed.u8(value(report.publication_status()));
  fixed.u8(value(report.rollback_status()));
  if (!fixed.valid() || fixed.size() != kReportSealEncodedSize)
    return 0U;
  return runtime::checkpoint_v2::crc64_ecma(fixed.data(), fixed.size());
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
void disseminate_failure_report(const runtime::MpiContext &mpi,
                                detail::CheckpointV2ReportValues &values,
                                CheckpointV2FailureReason reason,
                                CheckpointV2Phase phase,
                                int lowest_failing_rank,
                                std::uint64_t &collective_count,
                                std::string_view operation) {
  constexpr std::size_t count = 17U;
  std::array<std::uint64_t, count> common{};
  if (mpi.rank() == lowest_failing_rank) {
    std::uint64_t time_bits{};
    std::memcpy(&time_bits, &values.time_s, sizeof(time_bits));
    common = {static_cast<std::uint8_t>(reason),
              static_cast<std::uint8_t>(phase),
              values.step,
              time_bits,
              values.global_logical_bytes,
              values.global_actual_bytes,
              values.manifest_crc64,
              values.file_count,
              values.crc_check_count,
              static_cast<std::uint8_t>(values.manifest_crc),
              static_cast<std::uint8_t>(values.exact_size_eof),
              static_cast<std::uint8_t>(values.fingerprint),
              static_cast<std::uint8_t>(values.partition),
              static_cast<std::uint8_t>(values.transaction_entry),
              static_cast<std::uint8_t>(values.publication),
              static_cast<std::uint8_t>(values.rollback),
              static_cast<std::uint64_t>(lowest_failing_rank)};
  }
  runtime::check_mpi_result(
      MPI_Bcast(common.data(), static_cast<int>(common.size()), MPI_UINT64_T,
                lowest_failing_rank, mpi.comm()),
      operation);
  ++collective_count;
  values.disposition = CheckpointV2Disposition::failed;
  values.reason = static_cast<CheckpointV2FailureReason>(common[0U]);
  values.phase = static_cast<CheckpointV2Phase>(common[1U]);
  values.step = common[2U];
  std::memcpy(&values.time_s, &common[3U], sizeof(values.time_s));
  values.global_logical_bytes = common[4U];
  values.global_actual_bytes = common[5U];
  values.manifest_crc64 = common[6U];
  values.file_count = common[7U];
  values.crc_check_count = common[8U];
  values.manifest_crc = static_cast<CheckpointV2CheckStatus>(common[9U]);
  values.exact_size_eof = static_cast<CheckpointV2CheckStatus>(common[10U]);
  values.fingerprint = static_cast<CheckpointV2CheckStatus>(common[11U]);
  values.partition = static_cast<CheckpointV2CheckStatus>(common[12U]);
  values.transaction_entry = static_cast<CheckpointV2CheckStatus>(common[13U]);
  values.publication = static_cast<CheckpointV2CheckStatus>(common[14U]);
  values.rollback = static_cast<CheckpointV2CheckStatus>(common[15U]);
  values.lowest_failing_rank = static_cast<int>(common[16U]);
  values.collective_count = collective_count;
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
  return make_fingerprint(
      "hundun.checkpoint-v2.resolved-case.v1", [&](Encoder &encoder) {
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
        for (double item :
             {config.time.initial_dt_s, config.time.min_dt_s,
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
  return make_fingerprint(
      "hundun.checkpoint-v2.topology-local.v1", [&](Encoder &encoder) {
    encoder.i32(rank);
    encoder.i32(rank_count);
    append(encoder, decomposition.process_grid());
    append(encoder, topology.global_extent());
    append(encoder, topology.owned_global_box());
    encoder.u64(topology.global_cell_count());
    encoder.u64(topology.global_face_count());
        for (const auto axis :
             {mesh::FaceAxis::x, mesh::FaceAxis::y, mesh::FaceAxis::z})
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
          encoder.u8(topology.cell_ownership(cell) ==
                             mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
    }
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      encoder.u64(topology.global_face_id(face));
      const auto logical = topology.logical_face(face);
          encoder.u8(logical.axis == mesh::FaceAxis::x   ? 0U
                     : logical.axis == mesh::FaceAxis::y ? 1U
                                                         : 2U);
      append(encoder, logical.coordinate);
          encoder.u8(topology.face_ownership(face) ==
                             mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
      encoder.u64(topology.global_cell_id(topology.owner(face)));
      const auto neighbour = topology.neighbour(face);
      encoder.boolean(neighbour.has_value());
      if (neighbour)
        encoder.u64(topology.global_cell_id(*neighbour));
      const auto patch = topology.patch_id(face);
      encoder.boolean(patch.has_value());
      if (patch)
            encoder.u64(*patch);
      const auto pair = topology.periodic_pair(face);
      encoder.boolean(pair.has_value());
      if (pair)
        encoder.u64(*pair);
    }
    for (const auto &patch : topology.patches()) {
          encoder.u64(patch.stable_id());
      encoder.string(std::string(patch.name()));
      encoder.u8(patch.pairing_kind() == mesh::PatchPairingKind::none ? 0U
                                                                      : 1U);
      encoder.boolean(patch.paired_patch_id().has_value());
      if (patch.paired_patch_id())
            encoder.u64(*patch.paired_patch_id());
      encoder.u64(patch.local_faces().size());
      for (const auto face : patch.local_faces())
        encoder.u64(face);
    }
  });
}

std::uint64_t local_geometry_fingerprint(const mesh::MeshTopology &topology,
                                         const mesh::MeshGeometry &geometry,
                                         const config::FlowCaseConfig &config,
                                         int rank) {
  return make_fingerprint(
      "hundun.checkpoint-v2.geometry-local.v1", [&](Encoder &encoder) {
    encoder.i32(rank);
    encoder.u8(geometry.mapping_kind() == mesh::MappingKind::uniform_box
                       ? 0U
                       : 1U);
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
          const auto owner_area =
              geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
          append(encoder, owner_area);
      append(encoder,
                 runtime::Real3{-owner_area.x, -owner_area.y, -owner_area.z});
      encoder.f64(geometry.face_area_m2(face));
      encoder.f64(geometry.face_skewness(face));
      encoder.f64(geometry.face_non_orthogonality_degrees(face));
    }
  });
}

std::uint64_t
boundary_fingerprint(const boundary::BoundaryRegistry &boundaries) {
  return make_fingerprint(
      "hundun.checkpoint-v2.boundary.v1", [&](Encoder &encoder) {
    encoder.u64(boundaries.scalar_count());
        for (std::size_t scalar = 0; scalar < boundaries.scalar_count();
             ++scalar)
      encoder.string(std::string(boundaries.scalar_name(scalar)));
    encoder.boolean(boundaries.open_domain());
    encoder.boolean(boundaries.velocity_inlet_patch_id().has_value());
    if (boundaries.velocity_inlet_patch_id())
          encoder.u64(*boundaries.velocity_inlet_patch_id());
    encoder.boolean(boundaries.pressure_outlet_patch_id().has_value());
    if (boundaries.pressure_outlet_patch_id())
          encoder.u64(*boundaries.pressure_outlet_patch_id());
    for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
      const auto &item = boundaries.patch(patch);
          encoder.u64(item.stable_id());
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
            encoder.u64(*item.paired_patch_id());
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
  return make_fingerprint(
      "hundun.checkpoint-v2.field-schema.v1", [&](Encoder &encoder) {
    const auto ids = ordered_fields(fields);
    encoder.u64(ids.size());
    for (std::size_t index = 0; index < ids.size(); ++index) {
      const auto id = ids[index];
      const auto &descriptor = registry.descriptor(id);
      const std::uint8_t role =
          index < 5U ? static_cast<std::uint8_t>(index) : 5U;
      encoder.u8(role);
          encoder.u64(index < 5U ? 0U : index - 5U);
          encoder.u64(id);
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

std::uint64_t
local_layout_fingerprint(const runtime::StructuredDecomposition &decomposition,
                         const mesh::MeshTopology &topology, int rank,
                         int rank_count) {
  return make_fingerprint(
      "hundun.checkpoint-v2.local-layout.v1", [&](Encoder &encoder) {
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
          encoder.u8(topology.cell_ownership(cell) ==
                             mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
      encoder.u64(topology.global_cell_id(cell));
    }
    encoder.u64(topology.local_face_count());
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
          encoder.u8(topology.face_ownership(face) ==
                             mesh::EntityOwnership::owned
                         ? 0U
                         : 1U);
      encoder.u64(topology.global_face_id(face));
    }
  });
}

std::array<std::uint64_t, 5>
fingerprints(const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
             const mesh::MeshTopology &topology,
             const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const config::FlowCaseConfig &config,
    const runtime::FieldRegistry &registry, const FlowFieldIds &fields,
    std::uint64_t &collective_count) {
  std::uint64_t local_topology{};
  std::uint64_t local_geometry{};
  std::uint64_t resolved{};
  std::uint64_t boundary{};
  std::uint64_t schema{};
  bool local_prepared = true;
  try {
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::local_topology);
    local_topology = local_topology_fingerprint(decomposition, topology,
                                                mpi.rank(), mpi.size());
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::local_geometry);
    local_geometry =
        local_geometry_fingerprint(topology, geometry, config, mpi.rank());
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::resolved_case);
    resolved = resolved_fingerprint(config);
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::boundary_registry);
    boundary = boundary_fingerprint(boundaries);
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::field_schema);
    schema = field_schema_fingerprint(registry, fields);
  } catch (...) {
    local_prepared = false;
  }
  auto preparation = runtime::checkpoint_v2::converge_phase(
      mpi, local_prepared, collective_count,
      "MPI_Allreduce(Checkpoint local fingerprint preparation)");
  if (!preparation.ok)
    throw runtime::checkpoint_v2::CollectivePreparationError(
        preparation.failing_rank,
        "Checkpoint v2 local fingerprint preparation failed");
  const auto topology_parts = runtime::checkpoint_v2::allgather_u64(
      mpi, &local_topology, 1U, collective_count,
      "MPI_Allgather(Checkpoint topology fingerprints)");
  const auto geometry_parts = runtime::checkpoint_v2::allgather_u64(
      mpi, &local_geometry, 1U, collective_count,
      "MPI_Allgather(Checkpoint geometry fingerprints)");
  std::uint64_t topology_common{};
  std::uint64_t geometry_common{};
  bool common_prepared = true;
  try {
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::topology_common);
    topology_common = make_fingerprint(
        "hundun.checkpoint-v2.topology-common.v1", [&](Encoder &encoder) {
          encoder.i32(mpi.size());
          append(encoder, decomposition.process_grid());
          append(encoder, topology.global_extent());
          encoder.u64(topology_parts.size());
          for (int rank = 0; rank < mpi.size(); ++rank) {
            encoder.i32(rank);
            encoder.u64(topology_parts[static_cast<std::size_t>(rank)]);
          }
        });
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::geometry_common);
    geometry_common = make_fingerprint(
        "hundun.checkpoint-v2.geometry-common.v1", [&](Encoder &encoder) {
          encoder.u8(geometry.mapping_kind() == mesh::MappingKind::uniform_box
                         ? 0U
                         : 1U);
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
  } catch (...) {
    common_prepared = false;
  }
  preparation = runtime::checkpoint_v2::converge_phase(
      mpi, common_prepared, collective_count,
      "MPI_Allreduce(Checkpoint common fingerprint preparation)");
  if (!preparation.ok)
    throw runtime::checkpoint_v2::CollectivePreparationError(
        preparation.failing_rank,
        "Checkpoint v2 common fingerprint preparation failed");
  return {resolved, topology_common, geometry_common, boundary, schema};
}

ByteVector
encode_global_payload(AcceptedStepMetadata metadata,
                      const TimeControlState &time,
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
      order == 0U ? MomentumTimeOrder::backward_euler : MomentumTimeOrder::bdf2;
  result.time.schema_version = decoder.u32();
  result.time.accepted_step = decoder.u64();
  result.time.proposed_next_dt_s = decoder.f64();
  result.time.last_accepted_dt_s = decoder.f64();
  const auto last_order = decoder.u8();
  if (last_order > 1U)
    throw runtime::Error("Checkpoint v2 controller order is invalid");
  result.time.last_accepted_order = last_order == 0U
                                        ? MomentumTimeOrder::backward_euler
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
  const auto bytes = static_cast<std::uint64_t>(values.size() * sizeof(double));
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
  encoder.u32(
      static_cast<std::uint32_t>(fields.transported_cell_fields.size()));
  encoder.u32(static_cast<std::uint32_t>(
      2U * (5U + fields.transported_cell_fields.size())));
  logical_bytes = 0U;
  const auto layer = [&](std::uint8_t layer_id, const FlowLayerValues &values) {
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
                  runtime::FunctionSpace::face_value, 3U, values.face_velocity,
                  logical_bytes);
    append_record(encoder, layer_id, 4U, 0U, fields.face_mass_flux,
                  runtime::FunctionSpace::face_value, 1U, values.face_mass_flux,
                  logical_bytes);
    for (std::size_t index = 0; index < fields.transported_cell_fields.size();
         ++index)
      append_record(encoder, layer_id, 5U, static_cast<std::uint32_t>(index),
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
      layout.cell_interior_extent.y <= 0 || layout.cell_interior_extent.z <= 0)
    throw runtime::Error("Checkpoint v2 cell extent is invalid");
  std::size_t cells = runtime::checkpoint_v2::checked_product(
      static_cast<std::size_t>(layout.cell_interior_extent.x),
      static_cast<std::size_t>(layout.cell_interior_extent.y));
  cells = runtime::checkpoint_v2::checked_product(
      cells, static_cast<std::size_t>(layout.cell_interior_extent.z));
  const auto transported = state.fields().transported_cell_fields.size();
  if (transported > static_cast<std::size_t>(
                        std::numeric_limits<std::uint32_t>::max() / 2U - 5U))
    throw runtime::Error("Checkpoint v2 transported-field count is invalid");
  const auto records = runtime::checkpoint_v2::checked_product(
      2U, runtime::checkpoint_v2::checked_sum_u64(5U, transported));
  std::uint64_t one_layer_values{};
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values, static_cast<std::uint64_t>(cells));
  one_layer_values = runtime::checkpoint_v2::checked_sum_u64(
      one_layer_values, runtime::checkpoint_v2::checked_product(cells, 3U));
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
  const auto headers = runtime::checkpoint_v2::checked_product(records, 32U);
  const auto payload = runtime::checkpoint_v2::checked_sum_u64(
      32U, runtime::checkpoint_v2::checked_sum_u64(headers, logical));
  return {cells, records, logical, payload};
}

std::pair<FlowLayerValues, FlowLayerValues>
decode_rank_payload(const ByteVector &bytes, const FlowState &state,
                    std::uint64_t &logical_bytes) {
  const auto shape = expected_rank_payload_shape(state);
  if (bytes.size() != runtime::checkpoint_v2::checked_size(shape.payload_bytes))
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
        !std::all_of(
            layer.transported_cell_fields.front().begin(),
                     layer.transported_cell_fields.front().end(),
            [](double item) { return item > 0.0 && std::isfinite(item); }))
      return false;
  }
  return true;
}

std::string rank_filename(int rank) {
  std::array<char, 32> buffer{};
  const int count =
      std::snprintf(buffer.data(), buffer.size(), "rank-%06d.v2.bin", rank);
  if (count <= 0 || static_cast<std::size_t>(count) >= buffer.size())
    throw runtime::Error("Checkpoint v2 rank filename is invalid");
  return std::string(buffer.data(), static_cast<std::size_t>(count));
}

using Convergence = runtime::checkpoint_v2::CollectiveResult;

std::uint64_t sum_small_u64(const runtime::MpiContext &mpi, std::uint64_t local,
                            std::uint64_t &collective_count,
                            std::string_view operation) {
  std::uint64_t result{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &result, 1, MPI_UINT64_T, MPI_SUM, mpi.comm()),
      operation);
  ++collective_count;
  return result;
}

void broadcast_root_u64(const runtime::MpiContext &mpi, std::uint64_t *items,
                        std::size_t count, std::uint64_t &collective_count,
                        std::string_view operation) {
  if (items == nullptr ||
      count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error("Checkpoint v2 broadcast buffer is invalid");
  runtime::check_mpi_result(
      MPI_Bcast(items, static_cast<int>(count), MPI_UINT64_T, 0, mpi.comm()),
      operation);
  ++collective_count;
}
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
  return runtime::checkpoint_v2::converge_phase(mpi, local_ok, collective_count,
                                                operation);
}

std::uint64_t common_fingerprint(const std::array<std::uint64_t, 5> &items,
                                 int ranks, runtime::Int3 grid,
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
  return runtime::checkpoint_v2::exact_directory_inventory(directory, expected);
}

CheckpointV2FailureReason file_failure_reason(
    runtime::checkpoint_v2::NumericFileFailure failure) noexcept {
  return failure == runtime::checkpoint_v2::NumericFileFailure::filesystem
             ? CheckpointV2FailureReason::filesystem
             : CheckpointV2FailureReason::file_integrity;
}

runtime::checkpoint_v2::RankWrapper decode_authenticated_rank_wrapper(
    const ByteVector &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::int32_t expected_rank,
    std::int32_t expected_rank_count, std::uint64_t expected_payload_size) {
  if (bytes.size() !=
          runtime::checkpoint_v2::checked_size(expected_actual_size) ||
      expected_actual_size !=
          runtime::checkpoint_v2::checked_sum_u64(32U, expected_payload_size) ||
      runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size()) !=
          expected_crc)
    throw runtime::Error("Checkpoint v2 rank wrapper authentication failed");
  auto wrapper =
      runtime::checkpoint_v2::decode_rank_wrapper(bytes, expected_payload_size);
  if (wrapper.rank != expected_rank ||
      wrapper.rank_count != expected_rank_count)
    throw runtime::Error("Checkpoint v2 rank wrapper identity is invalid");
  return wrapper;
}

std::uint64_t expected_manifest_actual_size(
    std::uint64_t global_payload_size, std::uint64_t rank_count) {
  const auto rank_records = runtime::checkpoint_v2::checked_product(
      runtime::checkpoint_v2::checked_size(rank_count), 82U);
  const auto header_and_payload =
      runtime::checkpoint_v2::checked_sum_u64(84U, global_payload_size);
  const auto result = runtime::checkpoint_v2::checked_sum_u64(
      header_and_payload, static_cast<std::uint64_t>(rank_records));
  static_cast<void>(runtime::checkpoint_v2::checked_size(result));
  return result;
}

runtime::checkpoint_v2::Manifest decode_authenticated_manifest(
    const ByteVector &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::uint64_t expected_rank_count,
    std::uint64_t expected_global_payload_size) {
  const auto exact_size = expected_manifest_actual_size(
      expected_global_payload_size, expected_rank_count);
  if (expected_rank_count > std::numeric_limits<std::uint32_t>::max() ||
      expected_actual_size != exact_size ||
      bytes.size() !=
          runtime::checkpoint_v2::checked_size(expected_actual_size) ||
      runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size()) !=
          expected_crc)
    throw runtime::Error("Checkpoint v2 manifest authentication failed");
  return runtime::checkpoint_v2::decode_manifest(
      bytes, static_cast<std::uint32_t>(expected_rank_count),
      expected_global_payload_size);
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
bool same(const runtime::checkpoint_v2::ManifestRankRecord &left,
          const runtime::checkpoint_v2::ManifestRankRecord &right) noexcept {
  return left.rank == right.rank &&
         same(left.owned_box_begin, right.owned_box_begin) &&
         same(left.owned_box_end, right.owned_box_end) &&
         left.filename == right.filename &&
         left.logical_byte_size == right.logical_byte_size &&
         left.actual_byte_size == right.actual_byte_size &&
         left.crc64 == right.crc64 &&
         left.local_layout_fingerprint == right.local_layout_fingerprint;
}

bool same(const runtime::checkpoint_v2::Manifest &left,
          const runtime::checkpoint_v2::Manifest &right) noexcept {
  if (left.rank_count != right.rank_count ||
      !same(left.process_grid, right.process_grid) ||
      left.fingerprints != right.fingerprints ||
      left.global_payload != right.global_payload ||
      left.ranks.size() != right.ranks.size())
    return false;
  for (std::size_t index = 0; index < left.ranks.size(); ++index)
    if (!same(left.ranks[index], right.ranks[index]))
      return false;
  return true;
}

runtime::checkpoint_v2::CompletedMarker decode_authenticated_completed_marker(
    const ByteVector &bytes, std::uint64_t expected_manifest_actual_size,
    std::uint64_t expected_manifest_crc64,
    std::uint64_t expected_common_fingerprint) {
  if (bytes.size() != 40U)
    throw runtime::Error("Checkpoint v2 completed marker size is invalid");
  const auto marker = runtime::checkpoint_v2::decode_completed_marker(bytes);
  if (marker.manifest_actual_size != expected_manifest_actual_size ||
      marker.manifest_crc64 != expected_manifest_crc64 ||
      marker.common_fingerprint != expected_common_fingerprint)
    throw runtime::Error(
        "Checkpoint v2 completed marker authentication failed");
  return marker;
}
#endif

} // namespace

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
void set_checkpoint_v2_preparation_fault(CheckpointV2PreparationPoint point,
                                         std::uint32_t calls_before) noexcept {
  checkpoint_preparation_fault = static_cast<CheckpointPreparationPoint>(point);
  checkpoint_preparation_fault_calls_before = calls_before;
}
std::vector<std::uint8_t> checkpoint_v2_encode_global_payload_for_test(
    AcceptedStepMetadata metadata, const TimeControlState &time,
    const std::optional<IdealGasClosureState> &closure) {
  return encode_global_payload(metadata, time, closure);
}
std::vector<std::uint8_t>
checkpoint_v2_encode_rank_payload_for_test(const FlowState &state,
                                           std::uint64_t &logical_bytes) {
  return encode_rank_payload(state, logical_bytes);
}
bool checkpoint_v2_authenticate_global_payload_for_test(
    const std::vector<std::uint8_t> &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_size) noexcept {
  try {
    if (runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size()) !=
        expected_crc)
      return false;
    static_cast<void>(decode_global_payload(bytes, expected_size));
    return true;
  } catch (...) {
    return false;
  }
}
bool checkpoint_v2_authenticate_rank_payload_for_test(
    const std::vector<std::uint8_t> &bytes, const FlowState &state,
    std::uint64_t expected_crc) noexcept {
  try {
    if (runtime::checkpoint_v2::crc64_ecma(bytes.data(), bytes.size()) !=
        expected_crc)
      return false;
    std::uint64_t logical_bytes{};
    static_cast<void>(decode_rank_payload(bytes, state, logical_bytes));
    return true;
  } catch (...) {
    return false;
  }
}
bool checkpoint_v2_authenticate_rank_wrapper_for_test(
    const std::vector<std::uint8_t> &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::int32_t expected_rank,
    std::int32_t expected_rank_count,
    std::uint64_t expected_payload_size) noexcept {
  try {
    static_cast<void>(decode_authenticated_rank_wrapper(
        bytes, expected_crc, expected_actual_size, expected_rank,
        expected_rank_count, expected_payload_size));
    return true;
  } catch (...) {
    return false;
  }
}
bool checkpoint_v2_authenticate_manifest_for_test(
    const std::vector<std::uint8_t> &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size,
    const runtime::checkpoint_v2::Manifest &expected) noexcept {
  try {
    const auto decoded = decode_authenticated_manifest(
        bytes, expected_crc, expected_actual_size, expected.rank_count,
        expected.global_payload.size());
    return same(decoded, expected);
  } catch (...) {
    return false;
  }
}
bool checkpoint_v2_authenticate_manifest_limits_for_test(
    const std::vector<std::uint8_t> &bytes, std::uint64_t expected_crc,
    std::uint64_t expected_actual_size, std::uint64_t expected_rank_count,
    std::uint64_t expected_global_payload_size) noexcept {
  try {
    static_cast<void>(decode_authenticated_manifest(
        bytes, expected_crc, expected_actual_size, expected_rank_count,
        expected_global_payload_size));
    return true;
  } catch (...) {
    return false;
  }
}
bool checkpoint_v2_authenticate_completed_marker_for_test(
    const std::vector<std::uint8_t> &bytes,
    std::uint64_t expected_manifest_actual_size,
    std::uint64_t expected_manifest_crc64,
    std::uint64_t expected_common_fingerprint) noexcept {
  try {
    static_cast<void>(decode_authenticated_completed_marker(
        bytes, expected_manifest_actual_size, expected_manifest_crc64,
        expected_common_fingerprint));
    return true;
  } catch (...) {
    return false;
  }
}
std::uint64_t checkpoint_v2_field_schema_fingerprint_for_test(
    const runtime::FieldRegistry &registry, const FlowFieldIds &fields) {
  return field_schema_fingerprint(registry, fields);
}
} // namespace test
#endif

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
bool CheckpointV2ReadResult::ideal_gas_closure_state_available()
    const noexcept {
  return restored_ && closure_.has_value();
}
const IdealGasClosureState &
CheckpointV2ReadResult::ideal_gas_closure_state() const {
  if (!ideal_gas_closure_state_available())
    throw runtime::Error("Checkpoint v2 has no ideal-gas closure state");
  return *closure_;
}

CheckpointV2Report
detail::CheckpointV2Access::failed(CheckpointV2Operation operation, int rank,
                                   CheckpointV2FailureReason reason,
    CheckpointV2Phase phase) noexcept {
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
detail::CheckpointV2Access::make(CheckpointV2ReportValues values) noexcept {
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
    const TimeControlState &time, std::optional<IdealGasClosureState> closure,
    const std::filesystem::path &directory) {
  static_cast<void>(boundaries);
  detail::CheckpointV2ReportValues values;
  values.operation = CheckpointV2Operation::write;
  values.rank = mpi.rank();
  std::uint64_t collectives = 0U;
  const auto fail = [&](CheckpointV2FailureReason reason,
                        CheckpointV2Phase phase, int rank) {
    disseminate_failure_report(mpi, values, reason, phase, rank, collectives,
                               "MPI_Bcast(Checkpoint failure report)");
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
  Convergence path_status;
  try {
    path_status = path_agrees(mpi, directory, collectives);
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return fail(CheckpointV2FailureReason::state, CheckpointV2Phase::preflight,
                error.failing_rank());
  }
  if (!path_status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, path_status.failing_rank);
  try {
    if (!detail::FlowStateCheckpointAccess::live(state)) {
      local_reason = CheckpointV2FailureReason::invalid_input;
      throw runtime::Error("Checkpoint v2 FlowState is not live");
    }
    metadata = state.metadata();
    values.step = metadata.step;
    values.time_s = metadata.time_s;
    if (detail::FlowStateCheckpointAccess::attempt_active(state)) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 FlowState has an active attempt");
    }
    if (!geometry.compatible(topology) ||
        !same(topology.global_extent(), decomposition.global_extent()) ||
        !same(topology.owned_global_box(), decomposition.owned_box()) ||
        topology.local_face_count() !=
            detail::FlowStateCheckpointAccess::layout(state).face_count ||
        !same(decomposition.local_extent(),
              detail::FlowStateCheckpointAccess::layout(state)
                  .cell_interior_extent) ||
        !same(config.mesh.cells, topology.global_extent())) {
      local_reason = CheckpointV2FailureReason::layout;
      throw runtime::Error("Checkpoint v2 layout is incompatible");
    }
    if (config.schema_version != 2 ||
        config.simulation_type !=
            config::SimulationType::variable_density_flow) {
      local_reason = CheckpointV2FailureReason::invalid_input;
      throw runtime::Error("Checkpoint v2 configuration is invalid");
    }
    if (!detail::TimeControlStateCodec::semantically_valid(
            config.time, config.density_model, metadata, time)) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 time-control state is invalid");
    }
    const bool ideal = config.density_model == config::DensityModel::ideal_gas;
    if (closure.has_value() != ideal) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 closure presence is invalid");
    }
    if (closure &&
        (!(closure->thermodynamic_pressure_pa > 0.0) ||
         !std::isfinite(closure->thermodynamic_pressure_pa) ||
         closure->revision == std::numeric_limits<std::uint64_t>::max() ||
         (closure->mode == IdealGasPressureMode::closed_dynamic) !=
             closure->target_mass_kg.has_value() ||
         (closure->target_mass_kg &&
          (!(*closure->target_mass_kg > 0.0) ||
           !std::isfinite(*closure->target_mass_kg))))) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 closure state is invalid");
    }
    history = state.snapshot(FlowLayer::history);
    committed = state.snapshot(FlowLayer::committed);
    if (!valid_layer(history, config) || !valid_layer(committed, config)) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 physical state is invalid");
    }
  } catch (...) {
    local_ok = false;
  }
  auto status = converge(mpi, local_ok, collectives,
                         "MPI_Allreduce(Checkpoint semantic preflight)");
  if (!status.ok) {
    values.fingerprint = CheckpointV2CheckStatus::not_checked;
    values.partition = local_reason == CheckpointV2FailureReason::layout
                           ? CheckpointV2CheckStatus::failed
                           : CheckpointV2CheckStatus::not_checked;
    return fail(local_reason, CheckpointV2Phase::preflight,
                status.failing_rank);
  }

  bool payload_ok = true;
  local_reason = CheckpointV2FailureReason::state;
  try {
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::local_layout);
    local_layout = local_layout_fingerprint(decomposition, topology, mpi.rank(),
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
  } catch (...) {
    payload_ok = false;
  }
  status = converge(mpi, payload_ok, collectives,
                    "MPI_Allreduce(Checkpoint payload prepare)");
  if (!status.ok) {
    values.fingerprint = CheckpointV2CheckStatus::not_checked;
    values.partition = CheckpointV2CheckStatus::passed;
    return fail(local_reason, CheckpointV2Phase::preflight,
                status.failing_rank);
  }

  const auto &registry = detail::FlowStateCheckpointAccess::registry(state);
  try {
    identity = fingerprints(mpi, decomposition, topology, geometry, boundaries,
                            config, registry, state.fields(), collectives);
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    values.fingerprint = CheckpointV2CheckStatus::not_checked;
    values.partition = CheckpointV2CheckStatus::passed;
    return fail(CheckpointV2FailureReason::state, CheckpointV2Phase::preflight,
                error.failing_rank());
  }
  Encoder common_authority;
  bool authority_prepared = true;
  try {
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::common_authority);
    for (const auto fingerprint : identity)
      common_authority.u64(fingerprint);
    common_authority.i32(mpi.size());
    append(common_authority, decomposition.process_grid());
    common_authority.u64(global_payload.size());
    common_authority.raw(global_payload.data(), global_payload.size());
  } catch (...) {
    authority_prepared = false;
  }
  status = converge(mpi, authority_prepared, collectives,
                    "MPI_Allreduce(Checkpoint common authority preparation)");
  if (!status.ok)
    return fail(CheckpointV2FailureReason::state, CheckpointV2Phase::preflight,
                status.failing_rank);
  try {
    status = runtime::checkpoint_v2::opaque_bytes_agreement(
        mpi, common_authority.bytes(), collectives,
        "MPI_Allreduce(Checkpoint common authority)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return fail(CheckpointV2FailureReason::state, CheckpointV2Phase::preflight,
                error.failing_rank());
  }
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
  auto rank_write_reason = CheckpointV2FailureReason::filesystem;
  bool rank_write_integrity_failed = false;
  try {
    const auto verified = runtime::checkpoint_v2::write_verified_temporary(
        rank_temp, rank_wrapper);
    rank_written = verified.actual_size == rank_wrapper.size() &&
                   verified.crc64 == values.local_crc64;
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    rank_write_reason = file_failure_reason(error.failure());
    rank_write_integrity_failed =
        error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    rank_written = false;
  } catch (...) {
    rank_written = false;
  }
  status = converge(mpi, rank_written, collectives,
                    "MPI_Allreduce(Checkpoint rank temporary)");
  if (!status.ok) {
    values.crc_check_count =
        sum_small_u64(mpi, rank_written ? 1U : 0U, collectives,
                      "MPI_Allreduce(Checkpoint verified rank files)");
    values.rank_crc = rank_written ? CheckpointV2CheckStatus::passed
                      : rank_write_integrity_failed
                          ? CheckpointV2CheckStatus::failed
                                   : CheckpointV2CheckStatus::not_checked;
    values.exact_size_eof = rank_write_integrity_failed
                                ? CheckpointV2CheckStatus::failed
                                : CheckpointV2CheckStatus::not_checked;
    return fail(rank_write_reason, CheckpointV2Phase::rank_temporary_file,
                status.failing_rank);
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
    values.file_count =
        sum_small_u64(mpi, rank_published ? 1U : 0U, collectives,
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
      values.local_logical_bytes,
      values.local_actual_bytes,
      values.local_crc64,
      local_layout};
  std::vector<std::uint64_t> gathered_records;
  std::uint64_t global_logical{};
  std::uint64_t rank_actual{};
  try {
    gathered_records = runtime::checkpoint_v2::allgather_u64(
        mpi, local_record.data(), local_record.size(), collectives,
        "MPI_Allgather(Checkpoint rank records)");
    global_logical = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, values.local_logical_bytes, collectives,
        "MPI_Allreduce(Checkpoint logical bytes)");
    rank_actual = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, values.local_actual_bytes, collectives,
        "MPI_Allreduce(Checkpoint actual bytes)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::rank_publish, error.failing_rank());
  }

  std::uint64_t manifest_size = 0U;
  std::uint64_t manifest_crc = 0U;
  std::uint64_t marker_size = 0U;
  bool manifest_ok = true;
  bool manifest_verified = false;
  auto manifest_write_reason = CheckpointV2FailureReason::filesystem;
  bool manifest_integrity_failed = false;
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
        manifest.ranks.push_back({signed32(0U),
             {signed32(1U), signed32(2U), signed32(3U)},
             {signed32(4U), signed32(5U), signed32(6U)},
                                  rank_filename(rank),
                                  gathered_records[offset + 7U],
                                  gathered_records[offset + 8U],
                                  gathered_records[offset + 9U],
             gathered_records[offset + 10U]});
      }
      const auto bytes = runtime::checkpoint_v2::encode_manifest(manifest);
      const auto expected_size = expected_manifest_actual_size(
          static_cast<std::uint64_t>(global_payload.size()),
          static_cast<std::uint64_t>(mpi.size()));
      if (bytes.size() != expected_size)
        throw runtime::Error("Checkpoint v2 manifest size is invalid");
      const auto temp = directory / "manifest.v2.bin.tmp";
      const auto verified =
          runtime::checkpoint_v2::write_verified_temporary(temp, bytes);
      manifest_size = verified.actual_size;
      manifest_crc = verified.crc64;
      manifest_verified = true;
      runtime::checkpoint_v2::publish_no_overwrite(temp, directory /
                                                             "manifest.v2.bin");
    } catch (const runtime::checkpoint_v2::NumericFileError &error) {
      manifest_write_reason = file_failure_reason(error.failure());
      manifest_integrity_failed =
          error.failure() ==
          runtime::checkpoint_v2::NumericFileFailure::integrity;
      manifest_ok = false;
    } catch (...) {
      manifest_ok = false;
    }
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest)");
  const bool manifest_verified_common =
      sum_small_u64(mpi, manifest_verified ? 1U : 0U, collectives,
                    "MPI_Allreduce(Checkpoint manifest verification)") == 1U;
  if (!status.ok) {
    values.manifest_crc =
        manifest_integrity_failed  ? CheckpointV2CheckStatus::failed
        : manifest_verified_common ? CheckpointV2CheckStatus::passed
                              : CheckpointV2CheckStatus::not_checked;
    values.exact_size_eof = manifest_integrity_failed
                                ? CheckpointV2CheckStatus::failed
                                : CheckpointV2CheckStatus::passed;
    if (manifest_verified_common)
      ++values.crc_check_count;
    values.publication = CheckpointV2CheckStatus::failed;
    return fail(manifest_write_reason, CheckpointV2Phase::manifest,
                status.failing_rank);
  }
  const std::array<std::uint64_t, 2> local_manifest{
      mpi.rank() == 0 ? manifest_size : 0U,
      mpi.rank() == 0 ? manifest_crc : 0U};
  auto manifest_authority = local_manifest;
  broadcast_root_u64(mpi, manifest_authority.data(), manifest_authority.size(),
                     collectives, "MPI_Bcast(Checkpoint manifest authority)");
  manifest_size = manifest_authority[0U];
  manifest_crc = manifest_authority[1U];
  values.manifest_crc64 = manifest_crc;
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 1U;

  marker_size = 40U;
  if (rank_actual > std::numeric_limits<std::uint64_t>::max() - manifest_size ||
      rank_actual + manifest_size >
          std::numeric_limits<std::uint64_t>::max() - marker_size) {
    values.publication = CheckpointV2CheckStatus::not_checked;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::completed_marker, 0);
  }
  values.global_logical_bytes = global_logical;
  values.global_actual_bytes = rank_actual + manifest_size + marker_size;

  bool marker_ok = true;
  bool marker_verified = false;
  auto marker_write_reason = CheckpointV2FailureReason::filesystem;
  bool marker_integrity_failed = false;
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
    } catch (const runtime::checkpoint_v2::NumericFileError &error) {
      marker_write_reason = file_failure_reason(error.failure());
      marker_integrity_failed =
          error.failure() ==
          runtime::checkpoint_v2::NumericFileFailure::integrity;
      marker_ok = false;
    } catch (...) {
      marker_ok = false;
    }
  }
  status =
      converge(mpi, marker_ok, collectives, "MPI_Allreduce(Checkpoint marker)");
  const bool marker_verified_common =
      sum_small_u64(mpi, marker_verified ? 1U : 0U, collectives,
                    "MPI_Allreduce(Checkpoint marker verification)") == 1U;
  if (!status.ok) {
    if (marker_verified_common)
      ++values.crc_check_count;
    values.publication = CheckpointV2CheckStatus::failed;
    values.exact_size_eof = marker_integrity_failed
                                ? CheckpointV2CheckStatus::failed
                                : CheckpointV2CheckStatus::passed;
    values.fingerprint = marker_integrity_failed
                             ? CheckpointV2CheckStatus::failed
                             : CheckpointV2CheckStatus::passed;
    return fail(marker_write_reason, CheckpointV2Phase::completed_marker,
                status.failing_rank);
  }
  ++values.crc_check_count;
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

CheckpointV2ReadResult
read_checkpoint_v2(const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
                   const mesh::MeshTopology &topology,
                   const mesh::MeshGeometry &geometry,
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
    if (entered)
      values.rollback = CheckpointV2CheckStatus::passed;
    disseminate_failure_report(mpi, values, reason, phase, rank, collectives,
                               "MPI_Bcast(Checkpoint read failure report)");
    return detail::CheckpointV2Access::make_read(
        detail::CheckpointV2Access::make(values), {}, {}, false);
  };

  bool local_ok = true;
  CheckpointV2FailureReason local_reason{
      CheckpointV2FailureReason::invalid_input};
  Convergence path_status;
  try {
    path_status = path_agrees(mpi, directory, collectives);
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    return fail(CheckpointV2FailureReason::state, CheckpointV2Phase::preflight,
                error.failing_rank());
  }
  if (!path_status.ok)
    return fail(CheckpointV2FailureReason::invalid_input,
                CheckpointV2Phase::preflight, path_status.failing_rank);
  try {
    if (!detail::FlowStateCheckpointAccess::live(state)) {
      local_reason = CheckpointV2FailureReason::invalid_input;
      throw runtime::Error("Checkpoint v2 FlowState is not live");
    }
    const auto metadata = state.metadata();
    values.step = metadata.step;
    values.time_s = metadata.time_s;
    if (detail::FlowStateCheckpointAccess::attempt_active(state) ||
        !detail::FlowStateCheckpointAccess::diagnostic_identity_can_advance(
            state)) {
      local_reason = CheckpointV2FailureReason::state;
      throw runtime::Error("Checkpoint v2 FlowState transaction is invalid");
    }
    if (!geometry.compatible(topology) ||
        !same(topology.global_extent(), decomposition.global_extent()) ||
        !same(topology.owned_global_box(), decomposition.owned_box()) ||
        topology.local_face_count() !=
            detail::FlowStateCheckpointAccess::layout(state).face_count ||
        !same(decomposition.local_extent(),
                    detail::FlowStateCheckpointAccess::layout(state)
                  .cell_interior_extent)) {
      local_reason = CheckpointV2FailureReason::layout;
      throw runtime::Error("Checkpoint v2 read layout is incompatible");
    }
    if (config.schema_version != 2 ||
        config.simulation_type !=
            config::SimulationType::variable_density_flow) {
      local_reason = CheckpointV2FailureReason::invalid_input;
      throw runtime::Error("Checkpoint v2 read configuration is invalid");
    }
  } catch (...) {
    local_ok = false;
  }
  auto status = converge(mpi, local_ok, collectives,
               "MPI_Allreduce(Checkpoint read preflight)");
  if (!status.ok)
    return fail(local_reason, CheckpointV2Phase::preflight,
                status.failing_rank);

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
  auto marker_failure_phase = CheckpointV2Phase::marker_read;
  try {
    const auto marker_bytes = runtime::checkpoint_v2::read_regular_file_exact(
        directory / "COMPLETED", 40U);
    marker_exact = true;
    marker = runtime::checkpoint_v2::decode_completed_marker(marker_bytes);
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    marker_reason = file_failure_reason(error.failure());
    marker_exact_failed = error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    marker_ok = false;
  } catch (const std::bad_alloc &) {
    marker_reason = CheckpointV2FailureReason::state;
    marker_failure_phase = CheckpointV2Phase::restore_prepare;
    marker_ok = false;
  } catch (...) {
    marker_ok = false;
  }
  status = converge(mpi, marker_ok, collectives,
                    "MPI_Allreduce(Checkpoint marker read)");
  if (!status.ok) {
    values.exact_size_eof = marker_exact ? CheckpointV2CheckStatus::passed
                     : marker_exact_failed
                           ? CheckpointV2CheckStatus::failed
                           : CheckpointV2CheckStatus::not_checked;
    return fail(marker_reason, marker_failure_phase, status.failing_rank);
  }
  values.file_count = 1U;

  bool manifest_ok = true;
  bool manifest_exact = false;
  bool manifest_exact_failed = false;
  bool manifest_crc_checked = false;
  bool manifest_crc_match = false;
  auto manifest_reason = CheckpointV2FailureReason::file_integrity;
  auto manifest_failure_phase = CheckpointV2Phase::manifest_read;
  try {
    const auto global_size = expected_global_payload_size(
        config.density_model, boundaries.open_domain());
    const auto expected_manifest_size = expected_manifest_actual_size(
        global_size, static_cast<std::uint64_t>(mpi.size()));
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
    manifest = decode_authenticated_manifest(
        manifest_bytes, marker.manifest_crc64, marker.manifest_actual_size,
        static_cast<std::uint32_t>(mpi.size()), global_size);
    global = decode_global_payload(manifest.global_payload, global_size);
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    manifest_reason = file_failure_reason(error.failure());
    manifest_exact_failed =
        error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    manifest_ok = false;
  } catch (const std::bad_alloc &) {
    manifest_reason = CheckpointV2FailureReason::state;
    manifest_failure_phase = CheckpointV2Phase::restore_prepare;
    manifest_ok = false;
  } catch (...) {
    manifest_ok = false;
  }
  status = converge(mpi, manifest_ok, collectives,
                    "MPI_Allreduce(Checkpoint manifest read)");
  if (!status.ok) {
    values.manifest_crc =
        !manifest_crc_checked ? CheckpointV2CheckStatus::not_checked
        : manifest_crc_match  ? CheckpointV2CheckStatus::passed
                                    : CheckpointV2CheckStatus::failed;
    values.exact_size_eof = manifest_exact ? CheckpointV2CheckStatus::passed
                       : manifest_exact_failed
                             ? CheckpointV2CheckStatus::failed
                             : CheckpointV2CheckStatus::not_checked;
    values.crc_check_count = manifest_crc_match ? 1U : 0U;
    return fail(manifest_reason, manifest_failure_phase, status.failing_rank);
  }
  values.manifest_crc = CheckpointV2CheckStatus::passed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  values.file_count = 2U;
  values.crc_check_count = 1U;

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
    return fail(inventory_reason, CheckpointV2Phase::manifest_read,
                status.failing_rank);

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
  int fingerprint_preparation_failure = -1;
  try {
    const auto &registry = detail::FlowStateCheckpointAccess::registry(state);
    const auto expected =
        fingerprints(mpi, decomposition, topology, geometry, boundaries, config,
                     registry, state.fields(), collectives);
    fingerprint_ok = manifest.fingerprints == expected &&
        marker.common_fingerprint ==
            common_fingerprint(expected, mpi.size(),
                               decomposition.process_grid(),
                               manifest.global_payload);
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    fingerprint_preparation_failure = error.failing_rank();
    fingerprint_ok = false;
  } catch (...) {
    fingerprint_ok = false;
  }
  if (fingerprint_preparation_failure >= 0)
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare,
                fingerprint_preparation_failure);
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
  auto rank_failure_phase = CheckpointV2Phase::rank_read;
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
    values.local_crc64 = runtime::checkpoint_v2::crc64_ecma(rank_bytes.data(),
                                                            rank_bytes.size());
    rank_crc_checked = true;
    if (values.local_crc64 != record.crc64)
      throw runtime::Error("Checkpoint v2 rank CRC is invalid");
    rank_crc_match = true;
    const auto wrapper = decode_authenticated_rank_wrapper(
        rank_bytes, record.crc64, expected_actual, mpi.rank(), mpi.size(),
        shape.payload_bytes);
    std::uint64_t decoded_logical{};
    layers = decode_rank_payload(wrapper.payload, state, decoded_logical);
    if (decoded_logical != record.logical_byte_size)
      throw runtime::Error("Checkpoint v2 rank logical size is invalid");
    values.local_logical_bytes = decoded_logical;
  } catch (const runtime::checkpoint_v2::NumericFileError &error) {
    rank_reason = file_failure_reason(error.failure());
    rank_exact_failed = error.failure() ==
        runtime::checkpoint_v2::NumericFileFailure::integrity;
    rank_ok = false;
  } catch (const std::bad_alloc &) {
    rank_reason = CheckpointV2FailureReason::state;
    rank_failure_phase = CheckpointV2Phase::restore_prepare;
    rank_ok = false;
  } catch (...) {
    rank_ok = false;
  }
  status = converge(mpi, rank_ok, collectives,
                    "MPI_Allreduce(Checkpoint rank read)");
  if (!status.ok) {
    const auto exact_count =
        sum_small_u64(mpi, rank_exact ? 1U : 0U, collectives,
                      "MPI_Allreduce(Checkpoint exact rank files)");
    const auto exact_failure_count =
        sum_small_u64(mpi, rank_exact_failed ? 1U : 0U, collectives,
                      "MPI_Allreduce(Checkpoint invalid rank sizes)");
    const auto crc_count =
        sum_small_u64(mpi, rank_crc_match ? 1U : 0U, collectives,
                      "MPI_Allreduce(Checkpoint matching rank CRCs)");
    const auto complete_count =
        sum_small_u64(mpi, rank_ok ? 1U : 0U, collectives,
                      "MPI_Allreduce(Checkpoint complete rank files)");
    values.file_count =
        runtime::checkpoint_v2::checked_sum_u64(2U, complete_count);
    values.crc_check_count =
        runtime::checkpoint_v2::checked_sum_u64(1U, crc_count);
    values.rank_crc = !rank_crc_checked ? CheckpointV2CheckStatus::not_checked
                      : rank_crc_match  ? CheckpointV2CheckStatus::passed
                                : CheckpointV2CheckStatus::failed;
    values.exact_size_eof =
        exact_count == static_cast<std::uint64_t>(mpi.size())
            ? CheckpointV2CheckStatus::passed
        : exact_failure_count != 0U ? CheckpointV2CheckStatus::failed
                  : CheckpointV2CheckStatus::not_checked;
    return fail(rank_reason, rank_failure_phase, status.failing_rank);
  }
  values.rank_crc = CheckpointV2CheckStatus::passed;
  values.file_count = static_cast<std::uint64_t>(mpi.size()) + 2U;
  values.crc_check_count = static_cast<std::uint64_t>(mpi.size()) + 1U;

  bool state_ok =
      valid_layer(layers.first, config) && valid_layer(layers.second, config);
  if (config.density_model == config::DensityModel::ideal_gas && global.closure)
    state_ok = detail::validate_ideal_gas_restore_state(
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
    replacement.emplace(detail::FlowStateCheckpointAccess::prepare_replacement(
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
  std::uint64_t rank_actual{};
  try {
    values.global_logical_bytes = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, values.local_logical_bytes, collectives,
        "MPI_Allreduce(Checkpoint restored logical bytes)");
    rank_actual = runtime::checkpoint_v2::allreduce_sum_u64(
        mpi, values.local_actual_bytes, collectives,
        "MPI_Allreduce(Checkpoint restored actual bytes)");
  } catch (const runtime::checkpoint_v2::CollectivePreparationError &error) {
    values.publication = CheckpointV2CheckStatus::not_checked;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare, error.failing_rank());
  }
  bool success_boundary_ready = true;
  try {
    if (marker.manifest_actual_size >
            std::numeric_limits<std::uint64_t>::max() - 40U ||
        rank_actual >
            std::numeric_limits<std::uint64_t>::max() -
                (marker.manifest_actual_size + 40U) ||
        collectives == std::numeric_limits<std::uint64_t>::max()) {
      success_boundary_ready = false;
    } else {
      values.global_actual_bytes =
          rank_actual + marker.manifest_actual_size + 40U;
      values.collective_count = collectives + 1U;
    }
    inject_checkpoint_preparation_fault(
        CheckpointPreparationPoint::final_success_boundary);
  } catch (...) {
    success_boundary_ready = false;
  }
  status = converge(mpi, success_boundary_ready, collectives,
                    "MPI_Allreduce(Checkpoint restore success boundary)");
  if (!status.ok) {
    values.publication = CheckpointV2CheckStatus::not_checked;
    return fail(CheckpointV2FailureReason::state,
                CheckpointV2Phase::restore_prepare, status.failing_rank);
  }
  auto completed_report = detail::CheckpointV2Access::make(values);
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
